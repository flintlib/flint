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

void ifft_negacyclic(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w,
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
          mpn_neg(*t1, *t1, limbs + 1);
          SWAP_PTRS(ii[i], *t1);

          fft_adjust(*t2, ii[n+i], n - (n+i)/2, limbs, w);
          mpn_neg(*t2, *t2, limbs + 1);
          SWAP_PTRS(ii[n+i], *t2);

          i++;

          ifft_butterfly(*t1, *t2, ii[i], ii[n+i], i, limbs, w);
          SWAP_PTRS(ii[i],   *t1);
          SWAP_PTRS(ii[n+i], *t2);

          fft_adjust_sqrt2(*t1, ii[i], 2*n-i, limbs, w, *temp);
          mpn_neg(*t1, *t1, limbs + 1);
          SWAP_PTRS(ii[i], *t1);

          fft_adjust_sqrt2(*t2, ii[n+i], n-i, limbs, w, *temp);
          mpn_neg(*t2, *t2, limbs + 1);
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
          mpn_neg(*t1, *t1, limbs + 1);
          SWAP_PTRS(ii[i], *t1);

          fft_adjust(*t2, ii[n+i], n-i, limbs, w/2);
          mpn_neg(*t2, *t2, limbs + 1);
          SWAP_PTRS(ii[n+i], *t2);
       }
   }
}

