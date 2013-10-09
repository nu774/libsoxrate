/* 
 * Copyright (c) 2008 robs@users.sourceforge.net 
 * Copyright (c) 2011 nu774
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

#include <stdarg.h>
#include "libsoxrate.h"
#include "soxint.h"
#include "threaded_module.h"
#include "version.h"

struct lsx_fir_state_tag {
    int nchannels;
    sox_bool use_threads;
    lsx_thread_state_t thread;
    dft_filter_priv_t dft[1];
};

static int flow_channel(per_thread_state_t *pth);

lsx_fir_t *lsx_fir_create(unsigned nchannels, double *coefs, unsigned ncoefs,
			  unsigned post_peak, int use_threads)
{
    int i;
    lsx_fir_t *state;
    size_t off;

    if (!nchannels || !coefs || !ncoefs || ncoefs >= 0x10000)
	return 0;
    off = sizeof(dft_filter_priv_t) * (nchannels - 1);
    if ((state = (lsx_fir_t *)lsx_calloc(1, sizeof(lsx_fir_t) + off)) == 0)
	return 0;
    state->nchannels = nchannels;
    state->use_threads = (sox_bool)use_threads;
    for (i = 0; i < nchannels; ++i) {
	init(&state->dft[i], coefs, ncoefs, post_peak);
    }
    return state;
}

int lsx_fir_close(lsx_fir_t *state)
{
    int n;
    lsx_term_threads(&state->thread);
    for (n = 0; n < state->nchannels; ++n)
	stop(&state->dft[n]);
    free(state);
    return 0;
}

int lsx_fir_start(lsx_fir_t *state)
{
    int i;
    if (state->nchannels == 1)
	state->use_threads = sox_false;
    for (i = 0; i < state->nchannels; ++i)
	start(&state->dft[i]);
    if (lsx_init_threads(&state->thread, state->use_threads,
			 state->nchannels, 0, flow_channel) < 0)
	return -1;
    for (i = 0; i < state->nchannels; ++i)
	state->thread.pth[i].ctx = &state->dft[i];
    return 0;
}

int lsx_fir_process(lsx_fir_t *state, const float *ibuf, float *obuf,
		    size_t *ilen, size_t *olen)
{
    return lsx_process_threaded_interleaved(&state->thread, ibuf, obuf,
					    ilen, olen);
}

int lsx_fir_process_double(lsx_fir_t *state, const double *ibuf, double *obuf,
			   size_t *ilen, size_t *olen)
{
    return lsx_process_threaded_interleaved_double(&state->thread, ibuf, obuf,
						   ilen, olen);
}

int lsx_fir_process_noninterleaved(lsx_fir_t *state, const float * const *ibuf,
				   float **obuf, size_t *ilen, size_t *olen,
				   size_t istride, size_t ostride)
{
    return lsx_process_threaded_noninterleaved(&state->thread, ibuf, obuf,
					       ilen, olen, istride, ostride);
}

int lsx_fir_process_noninterleaved_double(lsx_fir_t *state,
					  const double * const *ibuf,
					  double **obuf, size_t *ilen,
					  size_t *olen, size_t istride,
					  size_t ostride)
{
    return lsx_process_threaded_noninterleaved_double(&state->thread, ibuf,
						      obuf, ilen, olen,
						      istride, ostride);
}

static int flow_channel(per_thread_state_t *pth)
{
    dft_filter_priv_t *dft = (dft_filter_priv_t*)pth->ctx;
    if (!pth->ilen)
	return drain(dft, pth->obuf._, &pth->olen);
    else
	return flow(dft, pth->ibuf._, pth->obuf._, &pth->ilen, &pth->olen);
}
