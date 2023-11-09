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

void butterfly_lshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1,
                       mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y)
{
   mp_limb_t cy, cy1, cy2;

   if (x == 0)
   {
      if (y == 0)
         cy = fft_sumdiff(t + x, u + y, i1, i2, limbs + 1);
      else
      {
         cy = fft_sumdiff(t, u + y, i1, i2, limbs - y);
         u[limbs] = -(cy&1);
         cy1 = cy>>1;
         cy = fft_sumdiff(t + limbs - y, u, i2 + limbs - y, i1 + limbs - y, y);
         t[limbs] = cy>>1;
         mpn_add_1(t + limbs - y, t + limbs - y, y + 1, cy1);
         cy1 = -(cy&1) + (i2[limbs] - i1[limbs]);
         mpn_addmod_2expp1_1(u + y, limbs - y, cy1);
         cy1 = -(i1[limbs] + i2[limbs]);
         mpn_addmod_2expp1_1(t, limbs, cy1);
      }
   } else if (y == 0)
   {
      cy = fft_sumdiff(t + x, u, i1, i2, limbs - x);
      t[limbs] = cy>>1;
      cy1 = cy&1;
      cy = fft_sumdiff(t, u + limbs - x, i1 + limbs - x, i2 + limbs - x, x);
      cy2 = mpn_neg(t, t, x);
      u[limbs] = -(cy&1);
      mpn_sub_1(u + limbs - x, u + limbs - x, x + 1, cy1);
      cy1 = -(cy>>1) - cy2;
      cy1 -= (i1[limbs] + i2[limbs]);
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
      cy1 = (i2[limbs] - i1[limbs]);
      mpn_addmod_2expp1_1(u, limbs, cy1);
   } else if (x > y)
   {
      cy = fft_sumdiff(t + x, u + y, i1, i2, limbs - x);
      t[limbs] = cy>>1;
      cy1 = cy&1;
      cy = fft_sumdiff(t, u + y + limbs - x, i1 + limbs - x, i2 + limbs - x, x - y);
      cy2 = mpn_neg(t, t, x - y);
      u[limbs] = -(cy&1);
      mpn_sub_1(u + y + limbs - x, u + y + limbs - x, x - y + 1, cy1);
      cy1 = (cy>>1) + cy2;
      cy = fft_sumdiff(t + x - y, u, i2 + limbs - y, i1 + limbs - y, y);
      cy2 = mpn_neg(t + x - y, t + x - y, y);
      cy1 = -(cy>>1) - mpn_sub_1(t + x - y, t + x - y, y, cy1) - cy2;
      cy1 -= (i1[limbs] + i2[limbs]);
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
      cy1 = -(cy&1) + (i2[limbs] - i1[limbs]);
      mpn_addmod_2expp1_1(u + y, limbs - y, cy1);
   } else if (x < y)
   {
      cy = fft_sumdiff(t + x, u + y, i1, i2, limbs - y);
      u[limbs] = -(cy&1);
      cy1 = cy>>1;
      cy = fft_sumdiff(t + x + limbs - y, u, i2 + limbs - y, i1 + limbs - y, y - x);
      t[limbs] = cy>>1;
      mpn_add_1(t + x + limbs - y, t + x + limbs - y, y - x + 1, cy1);
      cy1 = cy&1;
      cy = fft_sumdiff(t, u + y - x, i2 + limbs - x, i1 + limbs - x, x);
      cy1 = -(cy&1) - mpn_sub_1(u + y - x, u + y - x, x, cy1);
      cy1 += (i2[limbs] - i1[limbs]);
      mpn_addmod_2expp1_1(u + y, limbs - y, cy1);
      cy2 = mpn_neg(t, t, x);
      cy1 = -(cy>>1) - (i1[limbs] + i2[limbs]) - cy2;
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
   } else /* x == y */
   {
      cy = fft_sumdiff(t + x, u + x, i1, i2, limbs - x);
      t[limbs] = cy>>1;
      u[limbs] = -(cy&1);
      cy = fft_sumdiff(t, u, i2 + limbs - x, i1 + limbs - x, x);
      cy2 = mpn_neg(t, t, x);
      cy1 = -(cy>>1) - (i1[limbs] + i2[limbs]) - cy2;
      mpn_addmod_2expp1_1(t + x, limbs - x, cy1);
      cy1 = -(cy&1) + i2[limbs] - i1[limbs];
      mpn_addmod_2expp1_1(u + x, limbs - x, cy1);
  }
}

