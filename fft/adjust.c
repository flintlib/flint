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
      
void fft_adjust(mp_limb_t * r, mp_limb_t * i1, mp_size_t i, mp_size_t limbs, mp_bitcnt_t w)
{
   mp_bitcnt_t b1;
   mp_limb_t cy;
   mp_size_t x;

   b1 = i*w;
   x  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;

   if (x)
   {
      flint_mpn_copyi(r + x, i1, limbs - x);
      r[limbs] = 0;
      cy = mpn_neg_n(r, i1 + limbs - x, x);
      mpn_addmod_2expp1_1(r + x, limbs - x, -i1[limbs]);
      mpn_sub_1(r + x, r + x, limbs - x + 1, cy); 
      mpn_mul_2expmod_2expp1(r, r, limbs, b1);
   } else
      mpn_mul_2expmod_2expp1(r, i1, limbs, b1);
}
