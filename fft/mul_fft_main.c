/* 

Copyright 2009, 2011 William Hart. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL William Hart OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of William Hart.

*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
#include "ulong_extras.h"
#include "fft_tuning.h"

static int fft_tuning_table[5][2] = FFT_TAB;

void flint_mpn_mul_fft_main(mp_limb_t * r1, mp_limb_t * i1, mp_size_t n1, 
                        mp_limb_t * i2, mp_size_t n2)
{
   mp_size_t off, depth = 6;
   mp_size_t w = 1;
   mp_size_t n = ((mp_size_t) 1 << depth);
   mp_bitcnt_t bits = (n*w - (depth+1))/2;

   mp_bitcnt_t bits1 = n1*FLINT_BITS;
   mp_bitcnt_t bits2 = n2*FLINT_BITS;

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

