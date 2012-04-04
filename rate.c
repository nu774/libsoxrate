/* Effect: change sample rate     Copyright (c) 2008 robs@users.sourceforge.net
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

/* Inspired by, and builds upon some of the ideas presented in:
 * `The Quest For The Perfect Resampler' by Laurent De Soras;
 * http://ldesoras.free.fr/doc/articles/resampler-en.pdf */

#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "libsoxrate.h"
#include "soxint.h"
#include "fft4g.h"
#include "fifo.h"

#define  raw_coef_t double
#define  sample_t   double
#define  coef(coef_p, interp_order, fir_len, phase_num, coef_interp_num, fir_coef_num) coef_p[(fir_len) * ((interp_order) + 1) * (phase_num) + ((interp_order) + 1) * (fir_coef_num) + (interp_order - coef_interp_num)]

static sample_t * prepare_coefs(raw_coef_t const * coefs, int num_coefs,
    int num_phases, int interp_order, int multiplier)
{
  int i, j, length = num_coefs * num_phases;
  sample_t * result = lsx_malloc(length * (interp_order + 1) * sizeof(*result));
  double fm1 = coefs[0], f1 = 0, f2 = 0;

  for (i = num_coefs - 1; i >= 0; --i)
    for (j = num_phases - 1; j >= 0; --j) {
      double f0 = fm1, b = 0, c = 0, d = 0; /* = 0 to kill compiler warning */
      int pos = i * num_phases + j - 1;
      fm1 = (pos > 0 ? coefs[pos - 1] : 0) * multiplier;
      switch (interp_order) {
        case 1: b = f1 - f0; break;
        case 2: b = f1 - (.5 * (f2+f0) - f1) - f0; c = .5 * (f2+f0) - f1; break;
        case 3: c=.5*(f1+fm1)-f0;d=(1/6.)*(f2-f1+fm1-f0-4*c);b=f1-f0-d-c; break;
        default: if (interp_order) assert(0);
      }
      #define coef_coef(x) \
        coef(result, interp_order, num_coefs, j, x, num_coefs - 1 - i)
      coef_coef(0) = f0;
      if (interp_order > 0) coef_coef(1) = b;
      if (interp_order > 1) coef_coef(2) = c;
      if (interp_order > 2) coef_coef(3) = d;
      #undef coef_coef
      f2 = f1, f1 = f0;
    }
  return result;
}

typedef struct {    /* Data that are shared between channels and filters */
  sample_t   * poly_fir_coefs;
  dft_filter_t dft_filter[2];
} rate_shared_t;

struct stage;
typedef void (* stage_fn_t)(struct stage * input, fifo_t * output);
typedef struct stage {
  rate_shared_t * shared;
  fifo_t     fifo;
  int        pre;              /* Number of past samples to store */
  int        pre_post;         /* pre + number of future samples to store */
  int        preload;          /* Number of zero samples to pre-load the fifo */
  int        which;            /* Which, if any, of the 2 dft filters to use */
  stage_fn_t fn;
                               /* For poly_fir & spline: */
  union {                      /* 32bit.32bit fixed point arithmetic */
    #if defined(WORDS_BIGENDIAN)
    struct {int32_t integer; uint32_t fraction;} parts;
    #else
    struct {uint32_t fraction; int32_t integer;} parts;
    #endif
    int64_t all;
    #define MULT32 (65536. * 65536.)
  } at, step;
  int        L, remL, remM;
  double     out_in_ratio;
  fft_cache_t *cache;
} stage_t;

#define stage_occupancy(s) max(0, (int)fifo_occupancy(&(s)->fifo) - (s)->pre_post)
#define stage_read_p(s) ((sample_t *)fifo_read_ptr(&(s)->fifo) + (s)->pre)

static void cubic_spline_fn(stage_t * p, fifo_t * output_fifo)
{
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sample_t const * input = stage_read_p(p);
  sample_t * output = fifo_reserve(output_fifo, max_num_out);

  for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
    sample_t const * s = input + p->at.parts.integer;
    sample_t x = p->at.parts.fraction * (1 / MULT32);
    sample_t b = .5*(s[1]+s[-1])-*s, a = (1/6.)*(s[2]-s[1]+s[-1]-*s-4*b);
    sample_t c = s[1]-*s-a-b;
    output[i] = ((a*x + b)*x + c)*x + *s;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  fifo_read(&p->fifo, p->at.parts.integer, NULL);
  p->at.parts.integer = 0;
}

static void dft_stage_fn(stage_t * p, fifo_t * output_fifo)
{
  sample_t * output, tmp;
  int i, j, num_in = max(0, fifo_occupancy(&p->fifo));
  rate_shared_t const * s = p->shared;
  dft_filter_t const * f = &s->dft_filter[p->which];
  int const overlap = f->num_taps - 1;

  while (p->remL + p->L * num_in >= f->dft_length) {
    div_t divd = div(f->dft_length - overlap - p->remL + p->L - 1, p->L);
    sample_t const * input = fifo_read_ptr(&p->fifo);
    fifo_read(&p->fifo, divd.quot, NULL);
    num_in -= divd.quot;

    output = fifo_reserve(output_fifo, f->dft_length);
    if (p->L == 2 || p->L == 4) { /* F-domain */
      int portion = f->dft_length / p->L;
      memcpy(output, input, (unsigned)portion * sizeof(*output));
      lsx_safe_rdft(portion, 1, output, p->cache);
      for (i = portion + 2; i < (portion << 1); i += 2)
        output[i] = output[(portion << 1) - i],
        output[i+1] = -output[(portion << 1) - i + 1];
      output[portion] = output[1];
      output[portion + 1] = 0;
      output[1] = output[0];
      for (portion <<= 1; i < f->dft_length; i += portion, portion <<= 1) {
        memcpy(output + i, output, portion * sizeof(*output));
        output[i + 1] = 0;
      }
    } else {
      if (p->L == 1)
        memcpy(output, input, f->dft_length * sizeof(*output));
      else {
        memset(output, 0, f->dft_length * sizeof(*output));
        for (j = 0, i = p->remL; i < f->dft_length; ++j, i += p->L)
          output[i] = input[j];
        p->remL = p->L - 1 - divd.rem;
      }
      lsx_safe_rdft(f->dft_length, 1, output, p->cache);
    }
    output[0] *= f->coefs[0];
    if (p->step.parts.integer) {
      output[1] *= f->coefs[1];
      for (i = 2; i < f->dft_length; i += 2) {
        tmp = output[i];
        output[i  ] = f->coefs[i  ] * tmp - f->coefs[i+1] * output[i+1];
        output[i+1] = f->coefs[i+1] * tmp + f->coefs[i  ] * output[i+1];
      }
      lsx_safe_rdft(f->dft_length, -1, output, p->cache);
      if (p->step.parts.integer != 1) {
        for (j = 0, i = p->remM; i < f->dft_length - overlap; ++j, i += p->step.parts.integer)
          output[j] = output[i];
        p->remM = i - (f->dft_length - overlap);
        fifo_trim_by(output_fifo, f->dft_length - j);
      }
      else fifo_trim_by(output_fifo, overlap);
    }
    else { /* F-domain */
      for (i = 2; i < (f->dft_length >> 1); i += 2) {
        tmp = output[i];
        output[i  ] = f->coefs[i  ] * tmp - f->coefs[i+1] * output[i+1];
        output[i+1] = f->coefs[i+1] * tmp + f->coefs[i  ] * output[i+1];
      }
      output[1] = f->coefs[i] * output[i] - f->coefs[i+1] * output[i+1];
      lsx_safe_rdft(f->dft_length >> 1, -1, output, p->cache);
      fifo_trim_by(output_fifo, (f->dft_length + overlap) >> 1);
    }
  }
}

static void setup_dft_stage(rate_shared_t * shared, int which, stage_t * stage, int L, int M)
{
  stage->fn = dft_stage_fn;
  stage->preload = shared->dft_filter[which].post_peak / L;
  stage->remL    = shared->dft_filter[which].post_peak % L;
  stage->L = L;
  stage->step.parts.integer = M;
  stage->which = which;
}

static void init_dft_filter(rate_shared_t * p, unsigned which, int num_taps,
    sample_t const h[], double Fp, double Fc, double Fn, double att,
    int multiplier, double phase, sox_bool allow_aliasing, fft_cache_t *cache)
{
  dft_filter_t * f = &p->dft_filter[which];
  int dft_length, i;

  if (f->num_taps)
    return;
  if (h) {
    dft_length = lsx_set_dft_length(num_taps);
    f->coefs = lsx_calloc(dft_length, sizeof(*f->coefs));
    for (i = 0; i < num_taps; ++i)
      f->coefs[(i + dft_length - num_taps + 1) & (dft_length - 1)]
          = h[abs(num_taps / 2 - i)] / dft_length * 2 * multiplier;
    f->post_peak = num_taps / 2;
  }
  else {
    double * h2 = lsx_design_lpf(Fp, 1., 2., allow_aliasing, att, &num_taps, 0, -1);

    if (phase != 50)
      lsx_fir_to_phase(&h2, &num_taps, &f->post_peak, phase, cache);
    else {
      if (Fn == 4 && ((num_taps - 1) & 4)) { /* preserve phase */
        double * h3 = lsx_calloc(num_taps + 4, sizeof(*h3));
        memcpy(h3 + 2, h2, num_taps * sizeof(*h3));
        free(h2);
        h2 = h3;
        num_taps += 4;
      }
      f->post_peak = num_taps / 2;
    }

    dft_length = lsx_set_dft_length(num_taps);
    f->coefs = lsx_calloc(dft_length, sizeof(*f->coefs));
    for (i = 0; i < num_taps; ++i)
      f->coefs[(i + dft_length - num_taps + 1) & (dft_length - 1)]
          = h2[i] / dft_length * 2 * multiplier;
    free(h2);
  }
  assert(num_taps & 1);
  f->num_taps = num_taps;
  f->dft_length = dft_length;
  lsx_safe_rdft(dft_length, 1, f->coefs, cache);
}

#include "rate_filters.h"

typedef struct {
  double     factor;
  uint64_t   samples_in, samples_out;
  int        input_stage_num, output_stage_num;
  stage_t    * stages;
  fft_cache_t cache;
} rate_t;

#define pre_stage p->stages[-1]
#define frac_stage p->stages[level]
#define post_stage p->stages[level + have_frac_stage]
#define have_frac_stage (realM * fracL != 1)

typedef enum {Default = -1, Quick, Low, Medium, High, Very} quality_t;

static void rate_init(rate_t * p, rate_shared_t * shared, double factor,
    quality_t quality, int interp_order, double phase, double bandwidth,
    sox_bool allow_aliasing)
{
  int i, preL = 1, preM = 1, level = 0, fracL = 1, postL = 1, postM = 1;
  sox_bool upsample = sox_false;
  double realM = factor;

  assert(factor > 0);
  p->factor = factor;

  if (quality < Quick || quality > Very)
    quality = High;

  if (quality != Quick) while (sox_true) {
    const int max_divisor = 2048;      /* Keep coef table size ~< 500kb */
    double epsilon;
    upsample = realM < 1;
    for (i = realM, level = 0; i >>= 1; ++level); /* log base 2 */
    realM /= 1 << (level + !upsample);
    epsilon = fabs((uint32_t)(realM * MULT32 + .5) / (realM * MULT32) - 1);
    for (i = 2; i <= max_divisor && fracL == 1; ++i) {
      double try_d = realM * i;
      int try = try_d + .5;
      if (fabs(try / try_d - 1) <= epsilon) { /* N.B. beware of long doubles */
        if (try == i)
          realM = 1, fracL = 2, level += !upsample, upsample = sox_false;
        else realM = try, fracL = i;
      }
    }
    if (upsample) {
      if (postL == 1 && (realM != 1 || fracL > 5) && fracL / realM > 4) {
        realM = realM * (postL = min((fracL / realM), 4)) / fracL, fracL = 1;
        continue;
      }
      else if ((realM == 2 && fracL == 3) || (realM == 3 && fracL == 4))
        preL = fracL, preM = realM, fracL = realM = 1;
      else if (fracL < 6 && realM == 1)
        preL = fracL, fracL = 1;
      else if (quality > Low) {
        preL = 2;
        if (fracL % preL)
          realM *= preL;
        else fracL /= preL;
      }
    }
    else {
      if (fracL > 2) {
        int L = fracL, M = realM;
        for (i = level + 1; i && !(L & 1); L >>= 1, --i);
        if (L < 3 && (M <<= i) < 7) {
          preL = L, preM = M, realM = fracL = 1, level = 0, upsample = sox_true;
          break;
        }
      }
      postM = 2;
      if (fracL == 2)
        --fracL, postM -= !level, level -= !!level;
    }
    break;
  }
  memset(&p->cache, 0, sizeof(p->cache));
  p->stages = (stage_t *)lsx_calloc((size_t)level + 4, sizeof(*p->stages)) + 1;
  for (i = -1; i <= level + 1; ++i) {
    p->stages[i].shared = shared;
    p->stages[i].cache = &p->cache;
  }
  p->output_stage_num = level;

  frac_stage.step.all = realM * MULT32 + .5;
  frac_stage.out_in_ratio = MULT32 * fracL / frac_stage.step.all;

  if (quality == Quick) {
    frac_stage.fn = cubic_spline_fn;
    frac_stage.pre_post = max(3, frac_stage.step.parts.integer);
    frac_stage.preload = frac_stage.pre = 1;
    ++p->output_stage_num;
  }
  else if (have_frac_stage) {
    int n = (4 - (quality == Low)) * upsample + range_limit(quality, Medium, Very) - Medium;
    poly_fir_t const * f = &poly_firs[n];
    poly_fir1_t const * f1;

    if (f->num_coefs & 1) {
      if (fracL != 1 && (fracL & 1))
        fracL <<= 1, realM *= 2, frac_stage.step.all <<= 1;
      frac_stage.at.all = fracL * .5 * MULT32 + .5;
    }
    frac_stage.L = fracL;

    if (interp_order < 0)
      interp_order = quality > High;
    interp_order = fracL == 1? 1 + interp_order : 0;
    f1 = &f->interp[interp_order];

    if (!frac_stage.shared->poly_fir_coefs) {
      int phases = fracL == 1? (1 << f1->phase_bits) : fracL;
      int num_taps = f->num_coefs * phases - 1;
      raw_coef_t * coefs = lsx_design_lpf(
          f->pass, f->stop, 1., sox_false, f->att, &num_taps, phases, -1.);
      assert(num_taps == f->num_coefs * phases - 1);
      frac_stage.shared->poly_fir_coefs =
          prepare_coefs(coefs, f->num_coefs, phases, interp_order, 1);
      free(coefs);
    }
    frac_stage.fn = f1->fn;
    frac_stage.pre_post = f->num_coefs - 1;
    frac_stage.pre = 0;
    frac_stage.preload = frac_stage.pre_post >> 1;
    ++p->output_stage_num;
  }
  if (quality == Low && !upsample) {  /* dft is slower here, so */
    post_stage.fn = half_sample_low;       /* use normal convolution */
    post_stage.pre_post = 2 * (array_length(half_fir_coefs_low) - 1);
    post_stage.preload = post_stage.pre = post_stage.pre_post >> 1;
    ++p->output_stage_num;
  }
  else if (quality != Quick) {
    typedef struct {double bw, a;} filter_t;
    static filter_t const filters[] = {
      {.724, 100}, {.931, 110}, {.931, 125}, {.931, 170}};
    filter_t const * f = &filters[quality - Low];
    double att = allow_aliasing? (34./33)* f->a : f->a; /* negate att degrade */
    double bw = bandwidth? 1 - (1 - bandwidth / 100) / LSX_TO_3dB : f->bw;
    double min = 1 - (allow_aliasing? LSX_MAX_TBW0A : LSX_MAX_TBW0) / 100;
    double pass = bw * fracL / realM / 2;
    assert((size_t)(quality - Low) < array_length(filters));

    if (preL * preM != 1) {
      init_dft_filter(shared, 0, 0, 0, bw, 1., (double)max(preL, preM), att, preL, phase, allow_aliasing, &p->cache);
      setup_dft_stage(shared, 0, &pre_stage, preL, preM == 2 && !allow_aliasing? 0 : preM);
      --p->input_stage_num;
    }
    else if (level && have_frac_stage && (1 - pass) / (1 - bw) > 2)
      init_dft_filter(shared, 0, 0, NULL, max(pass, min), 1., 2., att, 1, phase, allow_aliasing, &p->cache);

    if (postL * postM != 1) {
      init_dft_filter(shared, 1, 0, 0,
          bw * (upsample? factor * postL / postM : 1),
          1., (double)(upsample? postL : postM), att, postL, phase, allow_aliasing, &p->cache);
      setup_dft_stage(shared, 1, &post_stage, postL, postM == 2 && !allow_aliasing? 0 : postM);
      ++p->output_stage_num;
    }
  }
  for (i = p->input_stage_num; i <= p->output_stage_num; ++i) {
    stage_t * s = &p->stages[i];
    if (i >= 0 && i < level - have_frac_stage) {
      s->fn = half_sample_25;
      s->pre_post = 2 * (array_length(half_fir_coefs_25) - 1);
      s->preload = s->pre = s->pre_post >> 1;
    }
    else if (level && i == level - 1) {
      if (shared->dft_filter[0].num_taps)
        setup_dft_stage(shared, 0, s, 1, (int)allow_aliasing << 1);
      else *s = post_stage;
    }
    fifo_create(&s->fifo, (int)sizeof(sample_t));
    memset(fifo_reserve(&s->fifo, s->preload), 0, sizeof(sample_t)*s->preload);
  }
}

static void rate_process(rate_t * p)
{
  stage_t * stage = p->stages + p->input_stage_num;
  int i;

  for (i = p->input_stage_num; i < p->output_stage_num; ++i, ++stage)
    stage->fn(stage, &(stage+1)->fifo);
}

static sample_t * rate_input(rate_t * p, sample_t const * samples, size_t n)
{
  p->samples_in += n;
  return fifo_write(&p->stages[p->input_stage_num].fifo, (int)n, samples);
}

static sample_t const * rate_output(rate_t * p, sample_t * samples, size_t * n)
{
  fifo_t * fifo = &p->stages[p->output_stage_num].fifo;
  p->samples_out += *n = min(*n, (size_t)fifo_occupancy(fifo));
  return fifo_read(fifo, (int)*n, samples);
}

static void rate_flush(rate_t * p)
{
  fifo_t * fifo = &p->stages[p->output_stage_num].fifo;
  uint64_t samples_out = p->samples_in / p->factor + .5;
  int64_t remaining = samples_out - p->samples_out;
  sample_t * buff = lsx_calloc(1024, sizeof(*buff));

  if ((int)remaining > 0) {
    while ((size_t)fifo_occupancy(fifo) < remaining) {
      rate_input(p, buff, (size_t) 1024);
      rate_process(p);
    }
    fifo_trim_to(fifo, (int)remaining);
    p->samples_in = 0;
  }
  free(buff);
}

static void rate_close(rate_t * p)
{
  rate_shared_t * shared = p->stages[0].shared;
  int i;

  for (i = p->input_stage_num; i <= p->output_stage_num; ++i)
    fifo_delete(&p->stages[i].fifo);
  free(shared->dft_filter[0].coefs);
  if (shared->dft_filter[1].coefs != shared->dft_filter[0].coefs)
    free(shared->dft_filter[1].coefs);
  free(shared->poly_fir_coefs);
  memset(shared, 0, sizeof(*shared));
  free(p->stages - 1);
  lsx_clear_fft_cache(&p->cache);
}

#include "rate_module.h"
