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
#include "mpn_extras.h"

void mul_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1, 
                        mp_srcptr i2, mp_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w)
{
   mp_size_t n = (UWORD(1)<<depth);
   flint_bitcnt_t bits1 = (n*w - (depth+1))/2; 
   
   mp_size_t r_limbs = n1 + n2;
   mp_size_t limbs = (n*w)/FLINT_BITS;
   mp_size_t size = limbs + 1;

   mp_size_t j1 = (n1*FLINT_BITS - 1)/bits1 + 1;
   mp_size_t j2 = (n2*FLINT_BITS - 1)/bits1 + 1;
   
   mp_size_t i, j, trunc;

   mp_limb_t ** ii, ** jj, * t1, * t2, * s1, * tt, * ptr;
   mp_limb_t c;
   
   ii = flint_malloc((4*(n + n*size) + 5*size)*sizeof(mp_limb_t));
   for (i = 0, ptr = (mp_limb_t *) ii + 4*n; i < 4*n; i++, ptr += size) 
   {
      ii[i] = ptr;
   }
   t1 = ptr;
   t2 = t1 + size;
   s1 = t2 + size;
   tt = s1 + size;
   
   if (i1 != i2)
   {
      jj = flint_malloc(4*(n + n*size)*sizeof(mp_limb_t));
      for (i = 0, ptr = (mp_limb_t *) jj + 4*n; i < 4*n; i++, ptr += size) 
      {
         jj[i] = ptr;
      }
   } else
      jj = ii;
   
   trunc = j1 + j2 - 1;
   if (trunc <= 2*n) trunc = 2*n + 1; /* trunc must be greater than 2n */
   trunc = 2*((trunc + 1)/2); /* trunc must be divisible by 2 */

   j1 = fft_split_bits(ii, i1, n1, bits1, limbs);
   for (j = j1 ; j < 4*n; j++)
      flint_mpn_zero(ii[j], limbs + 1);
   
   fft_truncate_sqrt2(ii, n, w, &t1, &t2, &s1, trunc);
    
   if (i1 != i2)
   {
      j2 = fft_split_bits(jj, i2, n2, bits1, limbs);
      for (j = j2 ; j < 4*n; j++)
         flint_mpn_zero(jj[j], limbs + 1);
      fft_truncate_sqrt2(jj, n, w, &t1, &t2, &s1, trunc);      
   } else j2 = j1;

   for (j = 0; j < trunc; j++)
   {
      mpn_normmod_2expp1(ii[j], limbs);
      if (i1 != i2) mpn_normmod_2expp1(jj[j], limbs);
      c = 2*ii[j][limbs] + jj[j][limbs];
      ii[j][limbs] = flint_mpn_mulmod_2expp1_basecase(ii[j], ii[j], jj[j], c, n*w, tt);
   }

   ifft_truncate_sqrt2(ii, n, w, &t1, &t2, &s1, trunc);
   for (j = 0; j < trunc; j++)
   {
      mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
      mpn_normmod_2expp1(ii[j], limbs);
   }
   
   flint_mpn_zero(r1, r_limbs);
   fft_combine_bits(r1, ii, j1 + j2 - 1, bits1, limbs, r_limbs);
     
   flint_free(ii);
   if (i1 != i2) flint_free(jj);
}

