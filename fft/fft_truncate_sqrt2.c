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
      
void fft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, 
                    mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, 
                         mp_size_t limbs, mp_bitcnt_t w, mp_limb_t * temp)
{
   mp_bitcnt_t wn = limbs*FLINT_BITS;
   mp_limb_t cy = 0;
   mp_size_t j = i/2, k = w/2;
   mp_size_t y;
   mp_bitcnt_t b1;
   int negate = 0;

   b1 = j + wn/4 + i*k;
   if (b1 >= wn) 
   {
      negate = 1;
      b1 -= wn;
   }
   y  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;
 
   /* sumdiff and multiply by 2^{j + wn/4 + i*k} */
   butterfly_lshB(s, t, i1, i2, limbs, 0, y);
   mpn_mul_2expmod_2expp1(t, t, limbs, b1);
   
   /* multiply by 2^{wn/2} */
   y = limbs/2;
   
   flint_mpn_copyi(temp + y, t, limbs - y);
   temp[limbs] = 0;
   if (y) cy = mpn_neg_n(temp, t + limbs - y, y);
   mpn_addmod_2expp1_1(temp + y, limbs - y, -t[limbs]);
   mpn_sub_1(temp + y, temp + y, limbs - y + 1, cy); 
   
   /* shift by an additional half limb (rare) */
   if (limbs & 1) 
       mpn_mul_2expmod_2expp1(temp, temp, limbs, FLINT_BITS/2);

   /* subtract */
   if (negate)
       mpn_sub_n(t, t, temp, limbs + 1);
   else
       mpn_sub_n(t, temp, t, limbs + 1);
}

void fft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
       mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w*n)/GMP_LIMB_BITS;
   
    if ((w & 1) == 0)
    {
        fft_truncate(ii, 2*n, w/2, t1, t2, trunc);
        return;
    }
   
    for (i = 0; i < trunc - 2*n; i++) 
    {   
        fft_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, limbs, w);
   
        SWAP_PTRS(ii[i],     *t1);
        SWAP_PTRS(ii[i+2*n], *t2);
 
        i++;
      
        fft_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, limbs, w, *temp);

        SWAP_PTRS(ii[i],     *t1);
        SWAP_PTRS(ii[2*n+i], *t2);
    }

    for (i = trunc - 2*n; i < 2*n; i++)
    {
        fft_adjust(ii[i+2*n], ii[i], i/2, limbs, w); 
         
        i++;

        fft_adjust_sqrt2(ii[i+2*n], ii[i], i, limbs, w, *temp); 
    }
   
    fft_radix2(ii, n, w, t1, t2);
    fft_truncate1(ii + 2*n, n, w, t1, t2, trunc - 2*n);
}
