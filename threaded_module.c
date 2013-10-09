#if defined(_MSC_VER) || defined(__MINGW32__)
#include <process.h>
#endif
#include <stdlib.h>
#include <memory.h>
#include "threaded_module.h"

#if defined(_MSC_VER) || defined(__MINGW32__)
static unsigned __stdcall worker_thread(void *arg);
#else
static void * worker_thread(void *arg);
#endif

#ifdef _MSC_VER
#define inline __inline
#define alloca _alloca
#endif

#ifdef _WIN32
inline void thread_join(pthread_t th)
{
    WaitForSingleObject(th, INFINITE);
    CloseHandle(th);
}

inline event_t event_new() { return CreateEventW(0, 0, 0, 0); }
inline void event_fire(event_t ev) { SetEvent(ev); }
inline int event_wait(event_t ev)
{
    return WaitForSingleObject(ev, INFINITE) == WAIT_OBJECT_0;
}
inline void event_close(event_t ev) { CloseHandle(ev); }
#else
static inline void thread_join(pthread_t th) { pthread_join(th, 0); }

static event_t event_new()
{
    event_t ev = lsx_malloc(sizeof(struct event_tag));
    if (!ev) return 0;
    pthread_mutex_init(&ev->mutex, 0);
    pthread_cond_init(&ev->cond, 0);
    ev->triggered = sox_false;
    return ev;
}

static void event_fire(event_t ev)
{
    pthread_mutex_lock(&ev->mutex);
    ev->triggered = sox_true;
    pthread_cond_signal(&ev->cond);
    pthread_mutex_unlock(&ev->mutex);
}

static int event_wait(event_t ev)
{
     pthread_mutex_lock(&ev->mutex);
     while (!ev->triggered)
         pthread_cond_wait(&ev->cond, &ev->mutex);
     ev->triggered = sox_false;
     pthread_mutex_unlock(&ev->mutex);
     return 1;
}

static void event_close(event_t ev)
{
    pthread_mutex_destroy(&ev->mutex);
    pthread_cond_destroy(&ev->cond);
    lsx_free(ev);
}
#endif

static inline float quantize(double v)
{
    const float anti_denormal = 1.0e-30f;
    float x = v;
    x += anti_denormal;
    x -= anti_denormal;
    return x;
}

#ifdef _WIN32
typedef unsigned (__stdcall *thread_proc)(void *);
#else
typedef void * (*thread_proc)(void *);
#endif

static int start_workers(lsx_thread_state_t *state, thread_proc func)
{
    int i;
    for (i = 0; i < state->count; ++i) {
	per_thread_state_t *ts = &state->pth[i];
	if (!(ts->evpro = event_new()))
	    return -1;
	if (!(ts->evcon = event_new()))
	    return -1;
#ifdef _WIN32
	if (!(ts->ht = (HANDLE)_beginthreadex(0, 0, func, ts, 0, &ts->tid)))
#else
	if (pthread_create(&ts->ht, 0, func, ts))
#endif
	    return -1;
    }
    return 0;
}

static int fire_and_wait_workers(lsx_thread_state_t *state, sox_bool join)
{
    int i, n = 0;

    for (i = 0; i < state->count; ++i) {
	per_thread_state_t *ts = &state->pth[i];
	if (ts->ht) event_fire(ts->evpro);
    }
    for (i = 0; i < state->count; ++i) {
	per_thread_state_t *ts = &state->pth[i];
	if (ts->ht) {
	    if (join)
		thread_join(ts->ht);
	    else
		event_wait(ts->evcon);
	}
    }
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
	if (ts->evpro) event_close(ts->evpro);
	if (ts->evcon) event_close(ts->evcon);
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
    lsx_free(state->pth);
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
    float **ivec = alloca(sizeof(float*) * state->count);
    float **ovec = alloca(sizeof(float*) * state->count);
    for (n = 0; n < state->count; ++n) {
	ivec[n] = ibuf + n;
	ovec[n] = obuf + n;
    }
    return lsx_process_threaded_noninterleaved(state, ivec, ovec,
					       ilen, olen,
					       state->count, state->count);
}

static unsigned roundup(unsigned n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

int lsx_process_threaded_noninterleaved(lsx_thread_state_t *state,
					const float * const *ibuf,
					float **obuf,
					size_t *ilen, size_t *olen,
					size_t istride, size_t ostride)
{
    int n;
    size_t i;
    unsigned count = ilen ? *ilen : 0;

    if (count > state->pth[0].ibuf.size) {
        for (n = 0; n < state->count; ++n) {
            unsigned size = state->pth[n].ibuf.size = roundup(count);
            state->pth[n].ibuf._ = realloc(state->pth[n].ibuf._,
                                           sizeof(double) * size);
        }
    }
    if (*olen > state->pth[0].obuf.size) {
        for (n = 0; n < state->count; ++n) {
            unsigned size = state->pth[n].obuf.size = roundup(*olen);
            state->pth[n].obuf._ = realloc(state->pth[n].obuf._,
                                           sizeof(double) * size);
        }
    }
    for (n = 0; n < state->count; ++n) {
	state->pth[n].ilen = count;
	state->pth[n].olen = *olen;
    }
    for (n = 0; n < state->count; ++n)
	for (i = 0; i < count; ++i)
	    state->pth[n].ibuf._[i] = ibuf[n][i * istride];

    if (run_filter(state) < 0)
	return -1;

    for (n = 0; n < state->count; ++n)
	for (i = 0; i < state->pth[0].olen; ++i)
	    obuf[n][i * ostride] = quantize(state->pth[n].obuf._[i]);
    if (ilen && *ilen)
	*ilen = state->pth[0].ilen;
    *olen = state->pth[0].olen;
    return 0;
}

int lsx_process_threaded_interleaved_double(lsx_thread_state_t *state,
					    const double *ibuf, double *obuf,
					    size_t *ilen, size_t *olen)
{
    int n;
    double **ivec = alloca(sizeof(double*) * state->count);
    double **ovec = alloca(sizeof(double*) * state->count);
    for (n = 0; n < state->count; ++n) {
	ivec[n] = ibuf + n;
	ovec[n] = obuf + n;
    }
    return lsx_process_threaded_noninterleaved_double(state, ivec, ovec,
						      ilen, olen,
						      state->count,
						      state->count);
}

int lsx_process_threaded_noninterleaved_double(lsx_thread_state_t *state,
					       const double * const *ibuf,
					       double **obuf, size_t *ilen,
					       size_t *olen, size_t istride,
					       size_t ostride)
{
    int n;
    size_t i;
    unsigned count = ilen ? *ilen : 0;

    if (count > state->pth[0].ibuf.size) {
        for (n = 0; n < state->count; ++n) {
            unsigned size = state->pth[n].ibuf.size = roundup(count);
            state->pth[n].ibuf._ = realloc(state->pth[n].ibuf._,
                                           sizeof(double) * size);
        }
    }
    if (*olen > state->pth[0].obuf.size) {
        for (n = 0; n < state->count; ++n) {
            unsigned size = state->pth[n].obuf.size = roundup(*olen);
            state->pth[n].obuf._ = realloc(state->pth[n].obuf._,
                                           sizeof(double) * size);
        }
    }
    for (n = 0; n < state->count; ++n) {
	state->pth[n].ilen = count;
	state->pth[n].olen = min(*olen, state->pth[n].obuf.size);
    }
    for (n = 0; n < state->count; ++n)
	for (i = 0; i < count; ++i)
	    state->pth[n].ibuf._[i] = ibuf[n][i * istride];

    if (run_filter(state) < 0)
	return -1;

    for (n = 0; n < state->count; ++n)
	for (i = 0; i < state->pth[0].olen; ++i)
	    obuf[n][i * ostride] = state->pth[n].obuf._[i];
    if (ilen && *ilen)
	*ilen = state->pth[0].ilen;
    *olen = state->pth[0].olen;
    return 0;
}

#ifdef _WIN32
static unsigned __stdcall worker_thread(void *arg)
#else
static void * worker_thread(void *arg)
#endif
{
    per_thread_state_t *pth = arg;
    while (event_wait(pth->evpro) && !pth->done) {
	if (pth->flow(pth) < 0)
	    pth->error = sox_true;
	event_fire(pth->evcon);
    }
    event_fire(pth->evcon);
    return 0;
}
