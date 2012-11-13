#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <pthread.h>
#endif
#include "soxint.h"

enum { IO_BUFSIZE = 4096 };

#ifdef _WIN32
typedef HANDLE event_t;
typedef HANDLE pthread_t;
#else
struct event_tag {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    sox_bool triggered;
};
typedef struct event_tag *event_t;
#endif

typedef struct per_thread_state_tag {
    pthread_t ht;
    event_t evpro, evcon;
    unsigned tid;
    sox_bool done;
    sox_bool error;

    void *ctx;
    int (*flow)(struct per_thread_state_tag *);
    size_t ilen, olen;
    double ibuf[IO_BUFSIZE], obuf[IO_BUFSIZE];
} per_thread_state_t;

typedef struct lsx_thread_state_tag {
    sox_bool use_threads;
    int count;
    per_thread_state_t *pth;
} lsx_thread_state_t;

int lsx_init_threads(lsx_thread_state_t *state, sox_bool use_threads,
		     int count, void **ctx, int (*flow)(per_thread_state_t *));
int lsx_term_threads(lsx_thread_state_t *state);
int lsx_process_threaded_interleaved(lsx_thread_state_t *state,
				     const float *ibuf, float *obuf,
				     size_t *ilen, size_t *olen);
int lsx_process_threaded_noninterleaved(lsx_thread_state_t *state,
					const float * const *ibuf,
					float **obuf,
					size_t *ilen, size_t *olen,
					size_t istride, size_t ostride);
int lsx_process_threaded_interleaved_double(lsx_thread_state_t *state,
					    const double *ibuf, double *obuf,
					    size_t *ilen, size_t *olen);
int lsx_process_threaded_noninterleaved_double(lsx_thread_state_t *state,
					       const double * const *ibuf,
					       double **obuf, size_t *ilen,
					       size_t *olen, size_t istride,
					       size_t ostride);
