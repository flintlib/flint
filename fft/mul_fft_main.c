/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
#include "ulong_extras.h"
#include "fft_tuning.h"

static int fft_tuning_table[5][2] = FFT_TAB;

void flint_mpn_mul_fft_main(mp_ptr r1, mp_srcptr i1, mp_size_t n1, 
                        mp_srcptr i2, mp_size_t n2)
{
   mp_size_t off, depth = 6;
   mp_size_t w = 1;
   mp_size_t n = ((mp_size_t) 1 << depth);
   flint_bitcnt_t bits = (n*w - (depth+1))/2;

   flint_bitcnt_t bits1 = n1*FLINT_BITS;
   flint_bitcnt_t bits2 = n2*FLINT_BITS;

   mp_size_t j1 = (bits1 - 1)/bits + 1;
   mp_size_t j2 = (bits2 - 1)/bits + 1;

   FLINT_ASSERT(n1 > 0);
   FLINT_ASSERT(n2 > 0);
   FLINT_ASSERT(j1 + j2 - 1 > 2*n);

   while (j1 + j2 - 1 > 4*n) /* find initial n, w */
   {
      if (w == 1) w = 2;
      else 
      {
         depth++;
         w = 1;
         n *= 2;
      }

      bits = (n*w - (depth+1))/2;
      j1 = (bits1 - 1)/bits + 1;
      j2 = (bits2 - 1)/bits + 1;
   }
   
   if (depth < 11)
   {
      mp_size_t wadj = 1;
      
      off = fft_tuning_table[depth - 6][w - 1]; /* adjust n and w */
      depth -= off;
      n = ((mp_size_t) 1 << depth);
      w *= ((mp_size_t) 1 << (2*off));
      
      if (depth < 6) wadj = ((mp_size_t) 1 << (6 - depth));

      if (w > wadj)
      {
         do { /* see if a smaller w will work */
            w -= wadj;
            bits = (n*w - (depth+1))/2;
            j1 = (bits1 - 1)/bits + 1;
            j2 = (bits2 - 1)/bits + 1;
         } while (j1 + j2 - 1 <= 4*n && w > wadj);  
         w += wadj;
      }

      mul_truncate_sqrt2(r1, i1, n1, i2, n2, depth, w);
   } else
   {
      if (j1 + j2 - 1 <= 3*n)
      {
         depth--;
         w *= 3;
      }
      mul_mfa_truncate_sqrt2(r1, i1, n1, i2, n2, depth, w);
   }
}

