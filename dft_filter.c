/* Abstract effect: dft filter     Copyright (c) 2008 robs@users.sourceforge.net
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
#include "soxint.h"
#include "fft4g.h"
#include "fifo.h"

typedef struct {
  size_t     samples_in, samples_out;
  fifo_t     input_fifo, output_fifo;
  dft_filter_t   filter, * filter_ptr;
  fft_cache_t  fft_cache;
} dft_filter_priv_t;

typedef dft_filter_t filter_t;
typedef dft_filter_priv_t priv_t;

static void init(dft_filter_priv_t *p, double *h, int n, int post_peak)
{
  int i;
  dft_filter_t *f = &p->filter;
  p->filter_ptr = f;
  if (!f->num_taps) {
    f->num_taps = n;
    f->post_peak = post_peak;
    f->dft_length = lsx_set_dft_length(f->num_taps);
    f->coefs = (double*)lsx_calloc(f->dft_length, sizeof(*f->coefs));
    for (i = 0; i < f->num_taps; ++i)
      f->coefs[(i + f->dft_length - f->num_taps + 1) & (f->dft_length - 1)]
	  = h[i] / f->dft_length * 2;
    lsx_safe_rdft(f->dft_length, 1, f->coefs, &p->fft_cache);
  }
}

static int start(dft_filter_priv_t * p)
{
  fifo_create(&p->input_fifo, (int)sizeof(double));
  memset(fifo_reserve(&p->input_fifo, p->filter_ptr->post_peak),
	 0, sizeof(double) * p->filter_ptr->post_peak);
  fifo_create(&p->output_fifo, (int)sizeof(double));
  return 0;
}

static void filter(dft_filter_priv_t * p)
{
  int i, num_in = max(0, fifo_occupancy(&p->input_fifo));
  filter_t const * f = p->filter_ptr;
  int const overlap = f->num_taps - 1;
  double * output;

  while (num_in >= f->dft_length) {
    double const * input = (double *)fifo_read_ptr(&p->input_fifo);
    fifo_read(&p->input_fifo, f->dft_length - overlap, NULL);
    num_in -= f->dft_length - overlap;

    output = (double *)fifo_reserve(&p->output_fifo, f->dft_length);
    fifo_trim_by(&p->output_fifo, overlap);
    memcpy(output, input, f->dft_length * sizeof(*output));

    lsx_safe_rdft(f->dft_length, 1, output, &p->fft_cache);
    output[0] *= f->coefs[0];
    output[1] *= f->coefs[1];
    for (i = 2; i < f->dft_length; i += 2) {
      double tmp = output[i];
      output[i  ] = f->coefs[i  ] * tmp - f->coefs[i+1] * output[i+1];
      output[i+1] = f->coefs[i+1] * tmp + f->coefs[i  ] * output[i+1];
    }
    lsx_safe_rdft(f->dft_length, -1, output, &p->fft_cache);
  }
}

static int flow(dft_filter_priv_t * p, const double * ibuf,
                double * obuf, size_t * isamp, size_t * osamp)
{
  size_t i, odone = min(*osamp, (size_t)fifo_occupancy(&p->output_fifo));
  double const * s = (double *)fifo_read(&p->output_fifo, (int)odone, NULL);

  for (i = 0; i < odone; ++i)
    *obuf++ = *s++;
  p->samples_out += odone;

  if (*isamp && odone < *osamp) {
    double * t = (double *)fifo_write(&p->input_fifo, (int)*isamp, NULL);
    p->samples_in += (int)*isamp;

    for (i = *isamp; i; --i)
      *t++ = *ibuf++;
    filter(p);
  }
  else *isamp = 0;
  *osamp = odone;
  return 0;
}

static int drain(dft_filter_priv_t * p, double * obuf, size_t * osamp)
{
  static size_t isamp = 0;
  size_t samples_out = p->samples_in;
  size_t remaining = samples_out - p->samples_out;
  double * buff = (double *)lsx_calloc(1024, sizeof(*buff));

  if ((int)remaining > 0) {
    while ((size_t)fifo_occupancy(&p->output_fifo) < remaining) {
      fifo_write(&p->input_fifo, 1024, buff);
      p->samples_in += 1024;
      filter(p);
    }
    fifo_trim_to(&p->output_fifo, (int)remaining);
    p->samples_in = 0;
  }
  free(buff);
  return flow(p, 0, obuf, &isamp, osamp);
}

static int stop(dft_filter_priv_t * p)
{
  fifo_delete(&p->input_fifo);
  fifo_delete(&p->output_fifo);
  free(p->filter_ptr->coefs);
  memset(p->filter_ptr, 0, sizeof(*p->filter_ptr));
  return 0;
}

#include "dft_module.h"
