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
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "fft4g.h"
#include "soxint.h"

static sox_bool is_power_of_2(int x)
{
  return !(x < 2 || (x & (x - 1)));
}

void lsx_clear_fft_cache(fft_cache_t *cp)
{
    cp->len = 0;
    free(cp->br);
    free(cp->sc);
    cp->br = 0;
    cp->sc = 0;
}

static void update_fft_cache(fft_cache_t *cp, int len)
{
  assert(is_power_of_2(len));
  assert(cp->len >= 0);
  if (len > cp->len) {
    int old_n = cp->len;
    cp->len = len;
    cp->br = lsx_realloc(cp->br, dft_br_len(len) * sizeof(cp->br[0]));
    cp->sc = lsx_realloc(cp->sc, dft_sc_len(len) * sizeof(cp->sc[0]));
    if (!old_n)
      cp->br[0] = 0;
  }
}

void lsx_safe_rdft(int len, int type, double * d, fft_cache_t *cp)
{
  update_fft_cache(cp, len);
  lsx_rdft(len, type, d, cp->br, cp->sc);
}

void lsx_safe_cdft(int len, int type, double * d, fft_cache_t *cp)
{
  update_fft_cache(cp, len);
  lsx_cdft(len, type, d, cp->br, cp->sc);
}
