/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fft.h"

void fft_adjust_sqrt2(mp_limb_t * r, mp_limb_t * i1,
            mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)
{
   flint_bitcnt_t wn = limbs*FLINT_BITS;
   mp_limb_t cy;
   mp_size_t j = i/2, k = w/2;
   mp_size_t y;
   flint_bitcnt_t b1;
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
      cy = mpn_neg(temp, i1 + limbs - y, y);
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
   if (y) cy = mpn_neg(temp, r + limbs - y, y);
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

