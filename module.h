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
#include <process.h>
#include <memory.h>
#define NOMINMAX
#include <windows.h>
#include "libsoxrate.h"
#include "version.h"

#define RATE_BUFSIZE 0x4000

typedef struct thread_state_tag {
    HANDLE ht;
    HANDLE evpro, evcon;
    unsigned tid;
    sox_bool done;
    sox_bool error;
} thread_state_t;

typedef struct per_channel_state_tag {
    thread_state_t th;
    rate_t rate;
    size_t ilen, olen;
    sample_t ibuf[RATE_BUFSIZE], obuf[RATE_BUFSIZE];
} per_channel_state_t;

struct lsx_rate_state_tag {
    int nchannels;
    sox_rate_t factor;
    quality_t quality;
    int phase;
    double bandwidth;
    sox_bool allow_aliasing;
    sox_bool use_threads;
    rate_shared_t shared;
    per_channel_state_t pcs[1];
};

static int start_workers(lsx_rate_t *state);
static int stop_workers(lsx_rate_t *state);
static int fire_and_wait_workers(lsx_rate_t *state, sox_bool join);
static unsigned __stdcall worker_thread(void *arg);

#define HANDLE_NO_MEMORY \
    GetExceptionCode() == STATUS_NO_MEMORY

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
    off = sizeof(per_channel_state_t) * (nchannels - 1);
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
    if (state->use_threads && stop_workers(state) < 0)
	return -1;
    for (n = 0; n < state->nchannels; ++n)
	rate_close(&state->pcs[n].rate);
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
    int i, n;
    if (state->nchannels == 1)
	state->use_threads = sox_false;
    __try {
	for (i = 0; i < state->nchannels; ++i)
	    rate_init(&state->pcs[i].rate, &state->shared, state->factor,
		      state->quality, -1, state->phase, state->bandwidth,
		      state->allow_aliasing);
	if (state->use_threads && start_workers(state) < 0) {
	    stop_workers(state);
	    for (n = 0; n < state->nchannels; ++n)
		rate_close(&state->pcs[n].rate);
	    return -1;
	}
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}

int lsx_rate_process(lsx_rate_t *state, const float *ibuf, float *obuf,
	size_t *ilen, size_t *olen)
{
    int n;
    size_t i, j = 0;

    __try {
	for (n = 0; n < state->nchannels; ++n) {
	    size_t count = ilen ? min(*ilen, RATE_BUFSIZE) : 0;
	    for (i = 0; i < count; ++i) 
		state->pcs[n].ibuf[i] = ibuf[state->nchannels * i + n];
	    state->pcs[n].ilen = count;
	    state->pcs[n].olen = min(*olen, RATE_BUFSIZE);
	    if (!count) rate_flush(&state->pcs[n].rate);
	}
	if (state->use_threads) {
	    if (fire_and_wait_workers(state, sox_false) < 0)
		return -1;
	    for (n = 0; n < state->nchannels; ++n)
		if (state->pcs[n].th.error)
		    return -1;
	} else {
	    for (n = 0; n < state->nchannels; ++n)
		if (flow_channel(&state->pcs[n]) < 0)
		    return -1;
	}
	for (i = 0; i < state->pcs[0].olen; ++i)
	    for (n = 0; n < state->nchannels; ++n)
		obuf[j++] = state->pcs[n].obuf[i];
	if (ilen && *ilen)
	    *ilen = state->pcs[0].ilen;
	*olen = state->pcs[0].olen;
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}

static int flow_channel(per_channel_state_t *pcs)
{
    __try {
	size_t odone = pcs->olen;
	rate_output(&pcs->rate, pcs->obuf, &odone);
	if (!pcs->ilen || odone == pcs->olen)
	    pcs->ilen = 0;
	else {
	    rate_input(&pcs->rate, pcs->ibuf, pcs->ilen);
	    rate_process(&pcs->rate);
	}
	pcs->olen = odone;
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}

static int start_workers(lsx_rate_t *state)
{
    int i;
    for (i = 0; i < state->nchannels; ++i) {
	thread_state_t *ts = &state->pcs[i].th;
	if (!(ts->evpro = CreateEventW(0, 0, 0, 0)))
	    return -1;
	if (!(ts->evcon = CreateEventW(0, 0, 0, 0)))
	    return -1;
	if (!(ts->ht = _beginthreadex(0, 0, worker_thread,
				&state->pcs[i], 0, &ts->tid)))
	    return -1;
    }
    return 0;
}

static int stop_workers(lsx_rate_t *state)
{
    int n;
    for (n = 0; n < state->nchannels; ++n)
	state->pcs[n].th.done = sox_true;
    if (fire_and_wait_workers(state, sox_true) < 0)
	return -1;
    for (n = 0; n < state->nchannels; ++n) {
	thread_state_t *ts = &state->pcs[n].th;
	if (ts->ht) CloseHandle(ts->ht);
	if (ts->evpro) CloseHandle(ts->evpro);
	if (ts->evcon) CloseHandle(ts->evcon);
	memset(ts, 0, sizeof(thread_state_t));
    }
    return 0;
}

static int fire_and_wait_workers(lsx_rate_t *state, sox_bool join)
{
    int i, n = 0;
    HANDLE *events;

    if ((events = _alloca(sizeof(HANDLE) * state->nchannels)) == 0)
	return -1;
    for (i = 0; i < state->nchannels; ++i) {
	thread_state_t *ts = &state->pcs[i].th;
	if (ts->ht) {
	    SetEvent(ts->evpro);
	    events[n++] = join ? ts->ht : ts->evcon;
	}
    }
    if (n) WaitForMultipleObjects(n, events, TRUE, INFINITE);
    return 0;
}

static unsigned __stdcall worker_thread(void *arg)
{
    per_channel_state_t *pcs = arg;
    while (WaitForSingleObject(pcs->th.evpro, INFINITE) == WAIT_OBJECT_0) {
	if (pcs->th.done)
	    break;
	if (flow_channel(pcs) < 0)
	    pcs->th.error = sox_true;
	SetEvent(pcs->th.evcon);
    }
    SetEvent(pcs->th.evcon);
    return 0;
}
