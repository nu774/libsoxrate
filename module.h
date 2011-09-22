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

#define RATE_BUFSIZE 4096

typedef struct thread_state_tag {
    HANDLE ht;
    HANDLE evpro, evcon;
    unsigned tid;
    sox_bool done;
} thread_state_t;

typedef struct per_channel_state_tag {
    thread_state_t th;
    rate_t rate;
    size_t ilen, olen;
    sample_t ibuf[RATE_BUFSIZE], obuf[RATE_BUFSIZE];
} per_channel_state_t;

struct lsx_rate_state_tag {
    quality_t quality;
    int phase;
    double bandwidth;
    sox_bool allow_aliasing;
    int nchannels;
    sox_rate_t factor;
    rate_shared_t shared;
    per_channel_state_t pcs[1];
};

static void start_workers(lsx_rate_t *state);
static void stop_workers(lsx_rate_t *state);
static void fire_and_wait_workers(lsx_rate_t *state);
static unsigned __stdcall worker_thread(void *arg);

lsx_rate_t *lsx_rate_create(unsigned nchannels,
	unsigned in_rate, unsigned out_rate)
{
    lsx_rate_t *state;
    size_t off;

    if (!nchannels || !in_rate || !out_rate)
	return 0;
    off = sizeof(per_channel_state_t) * (nchannels - 1);
    state = calloc(1, sizeof(lsx_rate_t) + off);
    if (!state)
	return 0;
    state->nchannels = nchannels;
    state->factor = (double)in_rate / out_rate;
    state->quality = Very;
    state->phase = 50; // linear;
    return state;
}

void lsx_rate_close(lsx_rate_t *state)
{
    int n;
    stop_workers(state);
    for (n = 0; n < state->nchannels; ++n)
	rate_close(&state->pcs[n].rate);
    free(state);
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
    } else if (type == SOX_RATE_ALLOW_ALIASING)
	state->allow_aliasing = !!va_arg(ap, int);
    else
	err = -1;
    va_end(ap);
    return err;
}

void lsx_rate_start(lsx_rate_t *state)
{
    int i;
    for (i = 0; i < state->nchannels; ++i)
	rate_init(&state->pcs[i].rate, &state->shared, state->factor,
		  state->quality, -1, state->phase, state->bandwidth,
		  state->allow_aliasing);
    start_workers(state);
}

size_t lsx_rate_process(lsx_rate_t *state, const float *ibuf, float *obuf,
	size_t *ilen, size_t *olen)
{
    int n;
    size_t i, j = 0;

    for (n = 0; n < state->nchannels; ++n) {
	size_t count = ilen ? min(*ilen, RATE_BUFSIZE) : 0;
	for (i = 0; i < count; ++i) 
	    state->pcs[n].ibuf[i] = ibuf[state->nchannels * i + n];
	state->pcs[n].ilen = count;
	state->pcs[n].olen = min(*olen, RATE_BUFSIZE);
	if (!count) rate_flush(&state->pcs[n].rate);
    }
    fire_and_wait_workers(state);
    for (i = 0; i < state->pcs[0].olen; ++i)
	for (n = 0; n < state->nchannels; ++n)
	    obuf[j++] = state->pcs[n].obuf[i];
    if (ilen && *ilen)
	*ilen = state->pcs[0].ilen;
    *olen = state->pcs[0].olen;
    return *olen;
}

static void flow_channel(per_channel_state_t *pcs)
{
    size_t odone = pcs->olen;
    rate_output(&pcs->rate, pcs->obuf, &odone);
    if (!pcs->ilen || odone == pcs->olen)
	pcs->ilen = 0;
    else {
	rate_input(&pcs->rate, pcs->ibuf, pcs->ilen);
	rate_process(&pcs->rate);
    }
    pcs->olen = odone;
}

static void start_workers(lsx_rate_t *state)
{
    int i;
    for (i = 0; i < state->nchannels; ++i) {
	thread_state_t *ts = &state->pcs[i].th;
	ts->evpro = CreateEventW(0, 0, 0, 0);
	ts->evcon = CreateEventW(0, 0, 0, 0);
	ts->ht = _beginthreadex(0, 0, worker_thread,
				&state->pcs[i], 0, &ts->tid);
    }
}

static void stop_workers(lsx_rate_t *state)
{
    int n;
    for (n = 0; n < state->nchannels; ++n)
	state->pcs[n].th.done = 1;
    fire_and_wait_workers(state);
    for (n = 0; n < state->nchannels; ++n) {
	thread_state_t *ts = &state->pcs[n].th;
	CloseHandle(ts->ht);
	CloseHandle(ts->evpro);
	CloseHandle(ts->evcon);
    }
}

static void fire_and_wait_workers(lsx_rate_t *state)
{
    int i;
    HANDLE *events;

    events = _alloca(sizeof(HANDLE) * state->nchannels);
    for (i = 0; i < state->nchannels; ++i) {
	SetEvent(state->pcs[i].th.evpro);
	events[i] = state->pcs[i].th.evcon;
    }
    WaitForMultipleObjects(state->nchannels, events, TRUE, INFINITE);
}

static unsigned __stdcall worker_thread(void *arg)
{
    per_channel_state_t *pcs = arg;
    while (WaitForSingleObject(pcs->th.evpro, INFINITE) == WAIT_OBJECT_0) {
	if (pcs->th.done)
	    break;
	flow_channel(pcs);
	SetEvent(pcs->th.evcon);
    }
    SetEvent(pcs->th.evcon);
    return 0;
}
