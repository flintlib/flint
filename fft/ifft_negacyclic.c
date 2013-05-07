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
      
void ifft_negacyclic(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w, 
                     mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/FLINT_BITS;

   ifft_radix2(ii, n/2, 2*w, t1, t2);
   ifft_radix2(ii+n, n/2, 2*w, t1, t2);

   if (w & 1)
   {
      for (i = 0; i < n; i++) 
      {   
          ifft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
          SWAP_PTRS(ii[i],   *t1);
          SWAP_PTRS(ii[n+i], *t2);

          fft_adjust(*t1, ii[i], n - i/2, limbs, w);
          mpn_neg_n(*t1, *t1, limbs + 1);
          SWAP_PTRS(ii[i], *t1);
            
          fft_adjust(*t2, ii[n+i], n - (n+i)/2, limbs, w);
          mpn_neg_n(*t2, *t2, limbs + 1);
          SWAP_PTRS(ii[n+i], *t2);
            
          i++;

          ifft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
          SWAP_PTRS(ii[i],   *t1);
          SWAP_PTRS(ii[n+i], *t2);

          fft_adjust_sqrt2(*t1, ii[i], 2*n-i, limbs, w, *temp);
          mpn_neg_n(*t1, *t1, limbs + 1);
          SWAP_PTRS(ii[i], *t1);
          
          fft_adjust_sqrt2(*t2, ii[n+i], n-i, limbs, w, *temp);
          mpn_neg_n(*t2, *t2, limbs + 1);
          SWAP_PTRS(ii[n+i], *t2);
       }
   } else
   {
       for (i = 0; i < n; i++) 
       {   
          ifft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
          SWAP_PTRS(ii[i],   *t1);
          SWAP_PTRS(ii[n+i], *t2);

          fft_adjust(*t1, ii[i], 2*n-i, limbs, w/2);
          mpn_neg_n(*t1, *t1, limbs + 1);
          SWAP_PTRS(ii[i], *t1);
            
          fft_adjust(*t2, ii[n+i], n-i, limbs, w/2);
          mpn_neg_n(*t2, *t2, limbs + 1);
          SWAP_PTRS(ii[n+i], *t2);
       }
   }
}

