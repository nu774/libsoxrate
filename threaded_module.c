#include <process.h>
#include <memory.h>
#include "threaded_module.h"

static unsigned __stdcall worker_thread(void *arg);

#define HANDLE_NO_MEMORY \
    GetExceptionCode() == STATUS_NO_MEMORY

static int start_workers(lsx_thread_state_t *state,
			 unsigned (__stdcall *func)(void *))
{
    int i;
    for (i = 0; i < state->count; ++i) {
	per_thread_state_t *ts = &state->pth[i];
	if (!(ts->evpro = CreateEventW(0, 0, 0, 0)))
	    return -1;
	if (!(ts->evcon = CreateEventW(0, 0, 0, 0)))
	    return -1;
	if (!(ts->ht = (HANDLE)_beginthreadex(0, 0, func, ts, 0, &ts->tid)))
	    return -1;
    }
    return 0;
}

static int fire_and_wait_workers(lsx_thread_state_t *state, sox_bool join)
{
    int i, n = 0;
    HANDLE *events;

    if ((events = _alloca(sizeof(HANDLE) * state->count)) == 0)
	return -1;
    for (i = 0; i < state->count; ++i) {
	per_thread_state_t *ts = &state->pth[i];
	if (ts->ht) {
	    SetEvent(ts->evpro);
	    events[n++] = join ? ts->ht : ts->evcon;
	}
    }
    if (n) WaitForMultipleObjects(n, events, TRUE, INFINITE);
    return 0;
}

static int stop_workers(lsx_thread_state_t *state)
{
    int n;
    for (n = 0; n < state->count ; ++n)
	state->pth[n].done = sox_true;
    if (fire_and_wait_workers(state, sox_true) < 0)
	return -1;
    for (n = 0; n < state->count; ++n) {
	per_thread_state_t *ts = &state->pth[n];
	if (ts->ht) CloseHandle(ts->ht);
	if (ts->evpro) CloseHandle(ts->evpro);
	if (ts->evcon) CloseHandle(ts->evcon);
	memset(ts, 0, sizeof(per_thread_state_t));
    }
    return 0;
}

int lsx_init_threads(lsx_thread_state_t *state, sox_bool use_threads,
		     int count, void **ctx, int (*flow)(per_thread_state_t *))
{
    int i;
    state->use_threads = use_threads;
    state->count = count;
    state->pth = lsx_calloc(count, sizeof(per_thread_state_t));
    for (i = 0; i < count; ++i)
	state->pth[i].flow = flow;
    if (ctx) {
	for (i = 0; i < count; ++i)
	    state->pth[i].ctx = ctx[i];
    }
    if (use_threads && start_workers(state, worker_thread) < 0) {
	stop_workers(state);
	return -1;
    }
    return 0;
}

int lsx_term_threads(lsx_thread_state_t *state)
{
    if (state->use_threads && stop_workers(state) < 0)
	return -1;
    free(state->pth);
    memset(state, 0, sizeof(lsx_thread_state_t));
    return 0;
}

static int run_filter(lsx_thread_state_t *state)
{
    int n;

    if (state->use_threads) {
	if (fire_and_wait_workers(state, sox_false) < 0)
	    return -1;
	for (n = 0; n < state->count; ++n)
	    if (state->pth[n].error)
		return -1;
    } else {
	for (n = 0; n < state->count; ++n)
	    if (state->pth[n].flow(&state->pth[n]) < 0)
		return -1;
    }
    return 0;
}

int lsx_process_threaded_interleaved(lsx_thread_state_t *state,
				     const float *ibuf, float *obuf,
				     size_t *ilen, size_t *olen)
{
    int n;
    size_t i, j = 0;

    __try {
	size_t count = ilen ? min(*ilen, IO_BUFSIZE) : 0;
	for (n = 0; n < state->count; ++n) {
	    state->pth[n].ilen = count;
	    state->pth[n].olen = min(*olen, IO_BUFSIZE);
	}
	j = 0;
	for (i = 0; i < count; ++i)
	    for (n = 0; n < state->count; ++n)
		state->pth[n].ibuf[i] = ibuf[j++];

	if (run_filter(state) < 0)
	    return -1;

	j = 0;
	for (i = 0; i < state->pth[0].olen; ++i)
	    for (n = 0; n < state->count; ++n)
		obuf[j++] = state->pth[n].obuf[i];
	if (ilen && *ilen)
	    *ilen = state->pth[0].ilen;
	*olen = state->pth[0].olen;
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}

int lsx_process_threaded_noninterleaved(lsx_thread_state_t *state,
					const float * const *ibuf,
					float **obuf,
					size_t *ilen, size_t *olen,
					size_t istride, size_t ostride)
{
    int n;
    size_t i;

    __try {
	size_t count = ilen ? min(*ilen, IO_BUFSIZE) : 0;
	for (n = 0; n < state->count; ++n) {
	    state->pth[n].ilen = count;
	    state->pth[n].olen = min(*olen, IO_BUFSIZE);
	}
	for (i = 0; i < count; ++i)
	    for (n = 0; n < state->count; ++n)
		state->pth[n].ibuf[i] = ibuf[n][i * istride];

	if (run_filter(state) < 0)
	    return -1;

	for (i = 0; i < state->pth[0].olen; ++i)
	    for (n = 0; n < state->count; ++n)
		obuf[n][i * ostride] = state->pth[n].obuf[i];
	if (ilen && *ilen)
	    *ilen = state->pth[0].ilen;
	*olen = state->pth[0].olen;
	return 0;
    } __except (HANDLE_NO_MEMORY) {
	return -1;
    }
}

static unsigned __stdcall worker_thread(void *arg)
{
    per_thread_state_t *pth = arg;
    while (WaitForSingleObject(pth->evpro, INFINITE) == WAIT_OBJECT_0) {
	if (pth->done)
	    break;
	if (pth->flow(pth) < 0)
	    pth->error = sox_true;
	SetEvent(pth->evcon);
    }
    SetEvent(pth->evcon);
    return 0;
}
