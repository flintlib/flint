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

void ifft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
   mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)
{
   flint_bitcnt_t wn = limbs*FLINT_BITS;
   mp_limb_t cy = 0;
   mp_size_t j = i/2, k = w/2;
   mp_size_t y2, y;
   mp_size_t b1;
   int negate = 1;

   b1 = wn - j - i*k - 1 + wn/4;
   if (b1 >= wn)
   {
      negate = 0;
      b1 -= wn;
   }
   y2 = b1/GMP_LIMB_BITS;
   b1 = b1%GMP_LIMB_BITS;

   /* multiply by small part of 2^{2*wn - j - ik - 1 + wn/4} */
   if (b1) mpn_mul_2expmod_2expp1(i2, i2, limbs, b1);

   /* multiply by 2^{wn/2} */
   y = limbs/2;

   flint_mpn_copyi(temp + y, i2, limbs - y);
   temp[limbs] = 0;
   if (y) cy = mpn_neg(temp, i2 + limbs - y, y);
   mpn_addmod_2expp1_1(temp + y, limbs - y, -i2[limbs]);
   mpn_sub_1(temp + y, temp + y, limbs - y + 1, cy);

   /* shift by an additional half limb (rare) */
   if (limbs & 1)
      mpn_mul_2expmod_2expp1(temp, temp, limbs, FLINT_BITS/2);

   /* subtract and negate... */
   if (negate) mpn_sub_n(i2, temp, i2, limbs + 1);
   else mpn_sub_n(i2, i2, temp, limbs + 1);

   /* ...negate and shift **left** by y2 limbs (i.e. shift right by
   (size - y2) limbs) and sumdiff */
   butterfly_rshB(s, t, i1, i2, limbs, 0, limbs - y2);
}

void ifft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
            mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc)
{
   mp_size_t i;
   mp_size_t limbs = (w*n)/FLINT_BITS;

   if ((w & 1) == 0)
   {
      ifft_truncate(ii, 2*n, w/2, t1, t2, trunc);
      return;
   }

   ifft_radix2(ii, n, w, t1, t2);

   for (i = trunc - 2*n; i < 2*n; i++)
   {
      fft_adjust(ii[i+2*n], ii[i], i/2, limbs, w);

      i++;

      fft_adjust_sqrt2(ii[i+2*n], ii[i], i, limbs, w, *temp);
   }

   ifft_truncate1(ii + 2*n, n, w, t1, t2, trunc - 2*n);

   for (i = 0; i < trunc - 2*n; i++)
   {
      ifft_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, limbs, w);

      SWAP_PTRS(ii[i], *t1);
      SWAP_PTRS(ii[2*n+i], *t2);

      i++;

      ifft_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, limbs, w, *temp);

      SWAP_PTRS(ii[i], *t1);
      SWAP_PTRS(ii[2*n+i], *t2);
   }

  for (i = trunc - 2*n; i < 2*n; i++)
     mpn_add_n(ii[i], ii[i], ii[i], limbs + 1);
}
