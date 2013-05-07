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
      
void fft_adjust_sqrt2(mp_limb_t * r, mp_limb_t * i1, 
            mp_size_t i, mp_size_t limbs, mp_bitcnt_t w, mp_limb_t * temp)
{
   mp_bitcnt_t wn = limbs*FLINT_BITS;
   mp_limb_t cy;
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
 
   /* multiply by 2^{j + wn/4 + i*k} */
   if (y)
   {
      flint_mpn_copyi(temp + y, i1, limbs - y);
      cy = mpn_neg_n(temp, i1 + limbs - y, y);
      temp[limbs] = 0;
      mpn_addmod_2expp1_1(temp + y, limbs - y, -i1[limbs]);
      mpn_sub_1(temp + y, temp + y, limbs - y + 1, cy); 
      mpn_mul_2expmod_2expp1(r, temp, limbs, b1);
   } else
      mpn_mul_2expmod_2expp1(r, i1, limbs, b1);

   /* multiply by 2^{wn/2} */
   y = limbs/2;
   cy = 0;

   flint_mpn_copyi(temp + y, r, limbs - y);
   temp[limbs] = 0;
   if (y) cy = mpn_neg_n(temp, r + limbs - y, y);
   mpn_addmod_2expp1_1(temp + y, limbs - y, -r[limbs]);
   mpn_sub_1(temp + y, temp + y, limbs - y + 1, cy); 
   
   /* shift by an additional half limb (rare) */
   if (limbs & 1) 
       mpn_mul_2expmod_2expp1(temp, temp, limbs, FLINT_BITS/2);

   /* subtract */
   if (negate)
      mpn_sub_n(r, r, temp, limbs + 1);
   else
      mpn_sub_n(r, temp, r, limbs + 1);
}

