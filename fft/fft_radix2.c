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
      
void fft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1, 
                   mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w)
{
   mp_size_t y;
   flint_bitcnt_t b1;

   b1 = i*w;
   y  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;
 
   butterfly_lshB(s, t, i1, i2, limbs, 0, y);
   mpn_mul_2expmod_2expp1(t, t, limbs, b1);
}

void fft_radix2(mp_limb_t ** ii, 
      mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2)
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
