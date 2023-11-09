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

void fft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t,
                    mp_limb_t * i1, mp_limb_t * i2, mp_size_t i,
                         mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)
{
   flint_bitcnt_t wn = limbs*FLINT_BITS;
   mp_limb_t cy = 0;
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

   /* sumdiff and multiply by 2^{j + wn/4 + i*k} */
   butterfly_lshB(s, t, i1, i2, limbs, 0, y);
   mpn_mul_2expmod_2expp1(t, t, limbs, b1);

   /* multiply by 2^{wn/2} */
   y = limbs/2;

   flint_mpn_copyi(temp + y, t, limbs - y);
   temp[limbs] = 0;
   if (y) cy = mpn_neg(temp, t + limbs - y, y);
   mpn_addmod_2expp1_1(temp + y, limbs - y, -t[limbs]);
   mpn_sub_1(temp + y, temp + y, limbs - y + 1, cy);

   /* shift by an additional half limb (rare) */
   if (limbs & 1)
       mpn_mul_2expmod_2expp1(temp, temp, limbs, FLINT_BITS/2);

   /* subtract */
   if (negate)
       mpn_sub_n(t, t, temp, limbs + 1);
   else
       mpn_sub_n(t, temp, t, limbs + 1);
}

void fft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
       mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc)
{
    mp_size_t i;
    mp_size_t limbs = (w*n)/GMP_LIMB_BITS;

    if ((w & 1) == 0)
    {
        fft_truncate(ii, 2*n, w/2, t1, t2, trunc);
        return;
    }

    for (i = 0; i < trunc - 2*n; i++)
    {
        fft_butterfly(*t1, *t2, ii[i], ii[2*n+i], i/2, limbs, w);

        SWAP_PTRS(ii[i],     *t1);
        SWAP_PTRS(ii[i+2*n], *t2);

        i++;

        fft_butterfly_sqrt2(*t1, *t2, ii[i], ii[2*n+i], i, limbs, w, *temp);

        SWAP_PTRS(ii[i],     *t1);
        SWAP_PTRS(ii[2*n+i], *t2);
    }

    for (i = trunc - 2*n; i < 2*n; i++)
    {
        fft_adjust(ii[i+2*n], ii[i], i/2, limbs, w);

        i++;

        fft_adjust_sqrt2(ii[i+2*n], ii[i], i, limbs, w, *temp);
    }

    fft_radix2(ii, n, w, t1, t2);
    fft_truncate1(ii + 2*n, n, w, t1, t2, trunc - 2*n);
}
