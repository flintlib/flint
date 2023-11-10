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

void fft_adjust(mp_limb_t * r, mp_limb_t * i1, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w)
{
   flint_bitcnt_t b1;
   mp_limb_t cy;
   mp_size_t x;

   b1 = i*w;
   x  = b1/FLINT_BITS;
   b1 = b1%FLINT_BITS;

   if (x)
   {
      flint_mpn_copyi(r + x, i1, limbs - x);
      r[limbs] = 0;
      cy = mpn_neg(r, i1 + limbs - x, x);
      mpn_addmod_2expp1_1(r + x, limbs - x, -i1[limbs]);
      mpn_sub_1(r + x, r + x, limbs - x + 1, cy);
      mpn_mul_2expmod_2expp1(r, r, limbs, b1);
   } else
      mpn_mul_2expmod_2expp1(r, i1, limbs, b1);
}
