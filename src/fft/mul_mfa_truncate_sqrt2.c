/* 
    Copyright (C) 2009, 2011, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fft.h"
#include "ulong_extras.h"

void mul_mfa_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
                        mp_srcptr i2, mp_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w)
{
   mp_size_t n = (UWORD(1)<<depth);
   flint_bitcnt_t bits1 = (n*w - (depth+1))/2; 
   mp_size_t sqrt = (UWORD(1)<<(depth/2));

   mp_size_t r_limbs = n1 + n2;
   mp_size_t limbs = (n*w)/FLINT_BITS;
   mp_size_t size = limbs + 1;

   mp_size_t j1 = (n1*FLINT_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*FLINT_BITS - 1)/bits1 + 1;
   
   mp_size_t i, j, trunc;

   mp_limb_t ** ii, ** jj, * ptr;
   mp_limb_t ** s1, ** t1, ** t2, ** tt;

   int N;

   TMP_INIT;

   TMP_START;

   N = flint_get_num_threads();
   ii = flint_malloc((4*(n + n*size) + 5*size*N)*sizeof(mp_limb_t));
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }

   s1 = TMP_ALLOC(N*sizeof(mp_limb_t *));
   t1 = TMP_ALLOC(N*sizeof(mp_limb_t *));
   t2 = TMP_ALLOC(N*sizeof(mp_limb_t *));
   tt = TMP_ALLOC(N*sizeof(mp_limb_t *));

   s1[0] = ptr;
   t1[0] = s1[0] + size*N;
   t2[0] = t1[0] + size*N;
   tt[0] = t2[0] + size*N;

   for (i = 1; i < N; i++)
   {
      s1[i] = s1[i - 1] + size;
      t1[i] = t1[i - 1] + size;
      t2[i] = t2[i - 1] + size;
      tt[i] = tt[i - 1] + 2*size;
   } 

   if (i1 != i2)
   {
      jj = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
      for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
      {
         jj[i] = ptr;
      }
   } else jj = ii;
   
   trunc = j1 + j2 - 1;
   if (trunc <= 2*n) trunc = 2*n + 1;
   trunc = 2*sqrt*((trunc + 2*sqrt - 1)/(2*sqrt)); /* trunc must be divisible by 2*sqrt */

   j1 = fft_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1 ; j < 4*n; j++)
      flint_mpn_zero(ii[j], limbs + 1);
   
   fft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);
   
   if (i1 != i2)
   {
      j2 = fft_split_bits(jj, i2, n2, bits1, limbs);
      for (j = j2 ; j < 4*n; j++)
         flint_mpn_zero(jj[j], limbs + 1);

      fft_mfa_truncate_sqrt2_outer(jj, n, w, t1, t2, s1, sqrt, trunc);
   } else j2 = j1;
   
   fft_mfa_truncate_sqrt2_inner(ii, jj, n, w, t1, t2, s1, sqrt, trunc, tt);
   ifft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);
       
   flint_mpn_zero(r1, r_limbs);
   fft_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   flint_free(ii);
   if (i1 != i2)
      flint_free(jj);

   TMP_END;
}
