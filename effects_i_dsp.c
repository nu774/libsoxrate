/* libSoX internal DSP functions.
 * All public functions & data are prefixed with lsx_ .
 *
 * Copyright (c) 2008 robs@users.sourceforge.net
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

#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <assert.h>
#include "soxint.h"
#include "fft4g.h"

static double lsx_bessel_I_0(double x)
{
  double term = 1, sum = 1, last_sum, x2 = x / 2;
  int i = 1;
  do {
    double y = x2 / i++;
    last_sum = sum, sum += term *= y * y;
  } while (sum != last_sum);
  return sum;
}

int lsx_set_dft_length(int num_taps) /* Set to 4 x nearest power of 2 */
{
  int result, n = num_taps;
  for (result = 8; n > 2; result <<= 1, n >>= 1);
  result = range_limit(result, 4096, 131072);
  assert(num_taps * 2 < result);
  return result;
}

static double lsx_kaiser_beta(double att)
{
  if (att > 100  ) return .1117 * att - 1.11;
  if (att > 50   ) return .1102 * (att - 8.7);
  if (att > 20.96) return .58417 * pow(att -20.96, .4) + .07886 * (att - 20.96);
  return 0;
}

static void lsx_apply_kaiser(double h[], const int num_points, double beta)
{
  int i, m = num_points - 1;
  for (i = 0; i <= m; ++i) {
    double x = 2. * i / m - 1;
    h[i] *= lsx_bessel_I_0(beta * sqrt(1 - x * x)) / lsx_bessel_I_0(beta);
  }
}

static double * lsx_make_lpf(int num_taps, double Fc, double beta, double scale, sox_bool dc_norm)
{
  int i, m = num_taps - 1;
  double * h = malloc(num_taps * sizeof(*h)), sum = 0;
  double mult = scale / lsx_bessel_I_0(beta);
  assert(Fc >= 0 && Fc <= 1);
  for (i = 0; i <= m / 2; ++i) {
    double x = M_PI * (i - .5 * m), y = 2. * i / m - 1;
    h[i] = x? sin(Fc * x) / x : Fc;
    sum += h[i] *= lsx_bessel_I_0(beta * sqrt(1 - y * y)) * mult;
    if (m - i != i)
      sum += h[m - i] = h[i];
  }
  for (i = 0; dc_norm && i < num_taps; ++i) h[i] *= scale / sum;
  return h;
}

static int lsx_lpf_num_taps(double att, double tr_bw, int k)
{                    /* TODO this could be cleaner, esp. for k != 0 */
  int n;
  if (att <= 80)
    n = .25 / M_PI * (att - 7.95) / (2.285 * tr_bw) + .5;
  else {
    double n160 = (.0425* att - 1.4) / tr_bw;   /* Half order for att = 160 */
    n = n160 * (16.556 / (att - 39.6) + .8625) + .5;  /* For att [80,160) */
  }
  return k? 2 * n : 2 * (n + (n & 1)) + 1; /* =1 %4 (0 phase 1/2 band) */
}

double * lsx_design_lpf(
    double Fp,      /* End of pass-band; ~= 0.01dB point */
    double Fc,      /* Start of stop-band */
    double Fn,      /* Nyquist freq; e.g. 0.5, 1, PI */
    sox_bool allow_aliasing,
    double att,     /* Stop-band attenuation in dB */
    int * num_taps, /* (Single phase.)  0: value will be estimated */
    int k)          /* Number of phases; 0 for single-phase */
{
  double tr_bw, beta;

  if (allow_aliasing)
    Fc += (Fc - Fp) * LSX_TO_3dB;
  Fp /= Fn, Fc /= Fn;        /* Normalise to Fn = 1 */
  tr_bw = LSX_TO_6dB * (Fc-Fp); /* Transition band-width: 6dB to stop points */

  if (!*num_taps)
    *num_taps = lsx_lpf_num_taps(att, tr_bw, k);
  beta = lsx_kaiser_beta(att);
  if (k)
    *num_taps = *num_taps * k - 1;
  else k = 1;
  return lsx_make_lpf(*num_taps, (Fc - tr_bw) / k, beta, (double)k, sox_true);
}

static double safe_log(double x)
{
  assert(x >= 0);
  if (x)
    return log(x);
  return -26;
}

void lsx_fir_to_phase(double * * h, int * len, int * post_len, double phase, fft_cache_t *cache)
{
  double * pi_wraps, * work, phase1 = (phase > 50 ? 100 - phase : phase) / 50;
  int i, work_len, begin, end, imp_peak = 0, peak = 0;
  double imp_sum = 0, peak_imp_sum = 0;
  double prev_angle2 = 0, cum_2pi = 0, prev_angle1 = 0, cum_1pi = 0;

  for (i = *len, work_len = 2 * 2 * 8; i > 1; work_len <<= 1, i >>= 1);

  work = calloc((size_t)work_len + 2, sizeof(*work)); /* +2: (UN)PACK */
  pi_wraps = malloc((((size_t)work_len + 2) / 2) * sizeof(*pi_wraps));

  memcpy(work, *h, *len * sizeof(*work));
  lsx_safe_rdft(work_len, 1, work, cache); /* Cepstral: */
  LSX_UNPACK(work, work_len);

  for (i = 0; i <= work_len; i += 2) {
    double angle = atan2(work[i + 1], work[i]);
    double detect = 2 * M_PI;
    double delta = angle - prev_angle2;
    double adjust = detect * ((delta < -detect * .7) - (delta > detect * .7));
    prev_angle2 = angle;
    cum_2pi += adjust;
    angle += cum_2pi;
    detect = M_PI;
    delta = angle - prev_angle1;
    adjust = detect * ((delta < -detect * .7) - (delta > detect * .7));
    prev_angle1 = angle;
    cum_1pi += fabs(adjust); /* fabs for when 2pi and 1pi have combined */
    pi_wraps[i >> 1] = cum_1pi;

    work[i] = safe_log(sqrt(sqr(work[i]) + sqr(work[i + 1])));
    work[i + 1] = 0;
  }
  LSX_PACK(work, work_len);
  lsx_safe_rdft(work_len, -1, work, cache);
  for (i = 0; i < work_len; ++i) work[i] *= 2. / work_len;

  for (i = 1; i < work_len / 2; ++i) { /* Window to reject acausal components */
    work[i] *= 2;
    work[i + work_len / 2] = 0;
  }
  lsx_safe_rdft(work_len, 1, work, cache);

  for (i = 2; i < work_len; i += 2) /* Interpolate between linear & min phase */
    work[i + 1] = phase1 * i / work_len * pi_wraps[work_len >> 1] +
        (1 - phase1) * (work[i + 1] + pi_wraps[i >> 1]) - pi_wraps[i >> 1];
        
  work[0] = exp(work[0]), work[1] = exp(work[1]);
  for (i = 2; i < work_len; i += 2) {
    double x = exp(work[i]);
    work[i    ] = x * cos(work[i + 1]);
    work[i + 1] = x * sin(work[i + 1]);
  }

  lsx_safe_rdft(work_len, -1, work, cache);
  for (i = 0; i < work_len; ++i) work[i] *= 2. / work_len;

  /* Find peak pos. */
  for (i = 0; i <= (int)(pi_wraps[work_len >> 1] / M_PI + .5); ++i) {
    imp_sum += work[i];          
    if (fabs(imp_sum) > fabs(peak_imp_sum)) {
      peak_imp_sum = imp_sum;
      peak = i;
    }
    if (work[i] > work[imp_peak]) /* For debug check only */
      imp_peak = i;
  }
  while (peak && fabs(work[peak-1]) > fabs(work[peak]) && work[peak-1] * work[peak] > 0)
    --peak;

  if (!phase1)
    begin = 0;
  else if (phase1 == 1)
    begin = peak - *len / 2;
  else {
    begin = (.997 - (2 - phase1) * .22) * *len + .5;
    end   = (.997 + (0 - phase1) * .22) * *len + .5;
    begin = peak - begin - (begin & 1);
    end   = peak + 1 + end + (end & 1);
    *len = end - begin;
    *h = realloc(*h, *len * sizeof(**h));
  }
  for (i = 0; i < *len; ++i) (*h)[i] =
    work[(begin + (phase > 50 ? *len - 1 - i : i) + work_len) & (work_len - 1)];
  *post_len = phase > 50 ? peak - begin : begin + *len - (peak + 1);

  free(pi_wraps), free(work);
}
