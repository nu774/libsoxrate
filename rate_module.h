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
#include "threaded_module.h"
#include "version.h"

struct lsx_rate_state_tag {
    int nchannels;
    sox_rate_t factor;
    quality_t quality;
    int phase;
    double bandwidth;
    sox_bool allow_aliasing;
    sox_bool use_threads;
    lsx_thread_state_t thread;
    rate_shared_t shared;
    rate_t rate[1];
};

#define HANDLE_NO_MEMORY \
    GetExceptionCode() == STATUS_NO_MEMORY

static int flow_channel(per_thread_state_t *pth);

const char *lsx_rate_version_string(void)
{
    return SOXRATE_VERSION_STRING;
}

lsx_rate_t *lsx_rate_create(unsigned nchannels,
	unsigned in_rate, unsigned out_rate)
{
    lsx_rate_t *state;
    size_t off;

    if (!nchannels || !in_rate || !out_rate)
	return 0;
    off = sizeof(rate_t) * (nchannels - 1);
    if ((state = lsx_calloc(1, sizeof(lsx_rate_t) + off)) == 0)
	return 0;
    state->nchannels = nchannels;
    state->factor = (double)in_rate / out_rate;
    state->quality = Very;
    state->phase = 50; // linear;
    state->use_threads = sox_true;
    return state;
}

int lsx_rate_close(lsx_rate_t *state)
{
    int n;
    lsx_term_threads(&state->thread);
    for (n = 0; n < state->nchannels; ++n)
	rate_close(&state->rate[n]);
    free(state);
    return 0;
}

int lsx_rate_config(lsx_rate_t *state, enum lsx_rate_config_e type, ...)
{
    int err = 0;
    va_list ap;

    va_start(ap, type);
    if (type == SOX_RATE_QUALITY) {
	int val = va_arg(ap, int);
	if (val < 0 || val > 2)
	    err = -1;
	else
	    state->quality = val + 2;
    } else if (type == SOX_RATE_PHASE_RESPONSE) {
	int val = va_arg(ap, int);
	if (val < 0 || val > 2)
	    err = -1;
	else
	    state->phase = val * 25;
    } else if (type == SOX_RATE_BANDWIDTH) {
	double val = va_arg(ap, double);
	if (val < 100 - LSX_MAX_TBW0A || val > 99.7)
	    err = -1;
	else
	    state->bandwidth = val;
    }
    else if (type == SOX_RATE_ALLOW_ALIASING)
	state->allow_aliasing = va_arg(ap, int);
    else if (type == SOX_RATE_USE_THREADS)
	state->use_threads = va_arg(ap, int);
    else
	err = -1;
    va_end(ap);
    return err;
}

int lsx_rate_start(lsx_rate_t *state)
{
    int i;
    if (state->nchannels == 1)
	state->use_threads = sox_false;
    __try {
	for (i = 0; i < state->nchannels; ++i)
	    rate_init(&state->rate[i], &state->shared, state->factor,
		      state->quality, -1, state->phase, state->bandwidth,
		      state->allow_aliasing);
	if (lsx_init_threads(&state->thread, state->use_threads,
			     state->nchannels, 0, flow_channel) < 0)
	    return -1;
	for (i = 0; i < state->nchannels; ++i)
	    state->thread.pth[i].ctx = &state->rate[i];
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}

int lsx_rate_process(lsx_rate_t *state, const float *ibuf, float *obuf,
	size_t *ilen, size_t *olen)
{
    return lsx_process_threaded_interleaved(&state->thread, ibuf, obuf,
					    ilen, olen);
}

static int flow_channel(per_thread_state_t *pth)
{
    __try {
	rate_t *rate = pth->ctx;
	size_t odone = pth->olen;
	if (!pth->ilen)
	    rate_flush(rate);
	rate_output(rate, pth->obuf, &odone);
	if (!pth->ilen || odone == pth->olen)
	    pth->ilen = 0;
	else {
	    rate_input(rate, pth->ibuf, pth->ilen);
	    rate_process(rate);
	}
	pth->olen = odone;
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}