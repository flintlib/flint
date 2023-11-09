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

void butterfly_rshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1,
                       mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y)
{
   mp_limb_t cy, cy1, cy2, cy3;

   if (x == 0)
   {
      if (y == 0)
      {
         cy = fft_sumdiff(t, u, i1, i2, limbs + 1);
      } else /* y != 0 */
      {
         cy = fft_sumdiff(t, u, i1, i2 + y, limbs - y);
         cy1 = (cy>>1);
         cy2 = -(cy&1);
         cy = fft_sumdiff(u + limbs - y, t + limbs - y, i1 + limbs - y, i2, y);
         u[limbs] = (cy>>1) + i1[limbs];
         t[limbs] = i1[limbs] - (cy&1);
         mpn_addmod_2expp1_1(t + limbs - y, y, cy1 + i2[limbs]);
         mpn_addmod_2expp1_1(u + limbs - y, y, cy2 - i2[limbs]);
      }
   } else if (y == 0) /* x != 0 */
   {
      cy = fft_sumdiff(t, u, i1 + x, i2, limbs - x);
      cy1 = (cy>>1);
      cy2 = -(cy&1);
      cy3 = mpn_neg(i1, i1, x);
      cy = fft_sumdiff(t + limbs - x, u + limbs - x, i1, i2 + limbs - x, x);
      u[limbs] = -cy3 - (cy&1) - i2[limbs];
      t[limbs] = -cy3 + i2[limbs] + (cy>>1);
      mpn_addmod_2expp1_1(t + limbs - x, x, cy1 + i1[limbs]);
      mpn_addmod_2expp1_1(u + limbs - x, x, cy2 + i1[limbs]);
   } else if (x == y)
   {
      cy = fft_sumdiff(t, u, i1 + x, i2 + x, limbs - x);
      cy1 = (cy>>1);
      cy2 = -(cy&1);
      cy = fft_sumdiff(t + limbs - x, u + limbs - x, i2, i1, x);
      cy3 = mpn_neg(t + limbs - x, t + limbs - x, x);
      u[limbs] = -(cy&1);
      t[limbs] = -(cy>>1) - cy3;
      mpn_addmod_2expp1_1(t + limbs - x, x, cy1 + i1[limbs] + i2[limbs]);
      mpn_addmod_2expp1_1(u + limbs - x, x, cy2 + i1[limbs] - i2[limbs]);
   } else if (x > y)
   {
      cy = fft_sumdiff(t + limbs - y, u + limbs - y, i2, i1 + x - y, y);
      cy3 = mpn_neg(t + limbs - y, t + limbs - y, y);
      t[limbs] = -(cy>>1) - cy3;
      u[limbs] = -(cy&1);
      cy3 = mpn_neg(i1, i1, x - y);
      cy = fft_sumdiff(t + limbs - x, u + limbs - x, i1, i2 + limbs - x + y, x - y);
      mpn_addmod_2expp1_1(t + limbs - y, y, (cy>>1) + i2[limbs] - cy3);
      mpn_addmod_2expp1_1(u + limbs - y, y, -(cy&1) - i2[limbs] - cy3);
      cy = fft_sumdiff(t, u, i1 + x, i2 + y, limbs - x);
      mpn_addmod_2expp1_1(t + limbs - x, x, (cy>>1) + i1[limbs]);
      mpn_addmod_2expp1_1(u + limbs - x, x, -(cy&1) + i1[limbs]);
   } else /* x < y */
   {
      cy = fft_sumdiff(t + limbs - x, u + limbs - x, i2 + y - x, i1, x);
      cy3 = mpn_neg(t + limbs - x, t + limbs - x, x);
      t[limbs] = -(cy>>1) - cy3;
      u[limbs] = -(cy&1);
      cy3 = mpn_neg(i2, i2, y - x);
      cy = fft_sumdiff(t + limbs - y, u + limbs - y, i1 + limbs - y + x, i2, y - x);
      mpn_addmod_2expp1_1(t + limbs - x, x, (cy>>1) + i1[limbs] - cy3);
      mpn_addmod_2expp1_1(u + limbs - x, x, -(cy&1) + i1[limbs] + cy3);
      cy = fft_sumdiff(t, u, i1 + x, i2 + y, limbs - y);
      mpn_addmod_2expp1_1(t + limbs - y, y, (cy>>1) + i2[limbs]);
      mpn_addmod_2expp1_1(u + limbs - y, y, -(cy&1) - i2[limbs]);
   }
}
