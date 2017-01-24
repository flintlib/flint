/* 
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "ulong_extras.h"
#include "fft.h"

void fft_mfa_truncate_sqrt2_inner(mp_limb_t ** ii, mp_limb_t ** jj, mp_size_t n, 
                   mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, 
                  mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc, mp_limb_t ** tt)
{
   mp_size_t i, j, s;
   mp_size_t n2 = (2*n)/n1;
   mp_size_t trunc2 = (trunc - 2*n)/n1;
   mp_size_t limbs = (n*w)/FLINT_BITS;
   mp_bitcnt_t depth = 0;
   mp_bitcnt_t depth2 = 0;
   int k = 0;

   while ((UWORD(1)<<depth) < n2) depth++;
   while ((UWORD(1)<<depth2) < n1) depth2++;

   ii += 2*n;
   jj += 2*n;

   /* convolutions on relevant rows */

#pragma omp parallel for private(s, i, j, k)
   for (s = 0; s < trunc2; s++)
   {
#if HAVE_OPENMP
      k = omp_get_thread_num();
#endif

      i = n_revbin(s, depth);
      fft_radix2(ii + i*n1, n1/2, w*n2, t1 + k, t2 + k);
      if (ii != jj) fft_radix2(jj + i*n1, n1/2, w*n2, t1 + k, t2 + k);
      
      for (j = 0; j < n1; j++)
      {
         mp_size_t t = i*n1 + j;
         mpn_normmod_2expp1(ii[t], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[t], limbs);
         fft_mulmod_2expp1(ii[t], ii[t], jj[t], n, w, tt[k]);
      }      
      
      ifft_radix2(ii + i*n1, n1/2, w*n2, t1 + k, t2 + k);
   }

   ii -= 2*n;
   jj -= 2*n;

   /* convolutions on rows */

#pragma omp parallel for private(i, j, k)
   for (i = 0; i < n2; i++)
   {
#if HAVE_OPENMP
      k = omp_get_thread_num();
#endif

      fft_radix2(ii + i*n1, n1/2, w*n2, t1 + k, t2 + k);
      if (ii != jj) fft_radix2(jj + i*n1, n1/2, w*n2, t1 + k, t2 + k);

      for (j = 0; j < n1; j++)
      {
         mp_size_t t = i*n1 + j;
         mpn_normmod_2expp1(ii[t], limbs);
         if (ii != jj) mpn_normmod_2expp1(jj[t], limbs);
         fft_mulmod_2expp1(ii[t], ii[t], jj[t], n, w, tt[k]);
      }      
      
      ifft_radix2(ii + i*n1, n1/2, w*n2, t1 + k, t2 + k);
   }
}

