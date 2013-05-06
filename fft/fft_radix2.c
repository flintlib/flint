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
      
void fft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1, 
                   mp_limb_t * i2, mp_size_t i, mp_size_t limbs, mp_bitcnt_t w)
{
   mp_size_t y;
   mp_bitcnt_t b1;

   b1 = i*w;
   y  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;
 
   butterfly_lshB(s, t, i1, i2, limbs, 0, y);
   mpn_mul_2expmod_2expp1(t, t, limbs, b1);
}

void fft_radix2(mp_limb_t ** ii, 
      mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/GMP_LIMB_BITS;
   
   if (n == 1) 
   {
      fft_butterfly(*t1, *t2, ii[0], ii[1], 0, limbs, w);

      SWAP_PTRS(ii[0], *t1);
      SWAP_PTRS(ii[1], *t2);
      
      return;
   }

   for (i = 0; i < n; i++) 
   {   
      fft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
   
      SWAP_PTRS(ii[i],   *t1);
      SWAP_PTRS(ii[n+i], *t2);
   }

   fft_radix2(ii,     n/2, 2*w, t1, t2);
   fft_radix2(ii + n, n/2, 2*w, t1, t2);
}
