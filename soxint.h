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
#ifndef SOXINT_H
#define SOXINT_H

#include <sys/types.h>

#define LSX_TO_6dB .5869
#define LSX_TO_3dB ((2/3.) * (.5 + LSX_TO_6dB))

#define LSX_MAX_TBW0 36.
#define LSX_MAX_TBW0A (LSX_MAX_TBW0 / (1 + LSX_TO_3dB))

#define sqr(a) ((a) * (a))

#define range_limit(x, lower, upper) (min(max(x, lower), upper))

/* Compile-time ("static") assertion */
/*   e.g. assert_static(sizeof(int) >= 4, int_type_too_small)    */
#define assert_static(e,f) enum {assert_static__##f = 1/(e)}
#define array_length(a) (sizeof(a)/sizeof(a[0]))

typedef enum sox_bool {
    sox_false,
    sox_true
} sox_bool;

typedef double sox_rate_t;

typedef struct fft_cache_tag {
    int len;
    int *br;
    double *sc;
} fft_cache_t;

typedef struct {
  int        dft_length, num_taps, post_peak;
  double     * coefs;
} dft_filter_t;

void lsx_clear_fft_cache(fft_cache_t *cp);
void lsx_safe_rdft(int len, int type, double * d, fft_cache_t *cp);
void lsx_safe_cdft(int len, int type, double * d, fft_cache_t *cp);
int lsx_set_dft_length(int num_taps);
double * lsx_design_lpf(
    double Fp,      /* End of pass-band; ~= 0.01dB point */
    double Fc,      /* Start of stop-band */
    double Fn,      /* Nyquist freq; e.g. 0.5, 1, PI */
    sox_bool allow_aliasing,
    double att,     /* Stop-band attenuation in dB */
    int * num_taps, /* (Single phase.)  0: value will be estimated */
    int k);         /* Number of phases; 0 for single-phase */
void lsx_fir_to_phase(double * * h, int * len, int * post_len, double phase0, fft_cache_t *cache);

void *lsx_malloc(size_t size);
void *lsx_calloc(size_t nelem, size_t size);
void *lsx_realloc(void *ptr, size_t size);

#endif
