/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

mp_size_t flint_mpn_gcd_full2(mp_ptr arrayg, 
                               mp_srcptr array1, mp_size_t limbs1,
		               mp_srcptr array2, mp_size_t limbs2, mp_ptr temp)
{
   mp_size_t s1 = 0, s2 = 0, m, b1, b2, mb, len1, len2, leng;
   mp_ptr in1, in2;
   mp_limb_t cy;

   /* find maximum power of 2 dividing inputs */
   b1 = mpn_scan1(array1 + s1, 0);
   b2 = mpn_scan1(array2 + s2, 0);
   
   /* get bit shifts [0, FLINT_BITS) and limb shifts */
   mb = FLINT_MIN(b1, b2) % FLINT_BITS;
   s1 = b1 / FLINT_BITS; b1 = b1 % FLINT_BITS; len1 = limbs1 - s1;
   s2 = b2 / FLINT_BITS; b2 = b2 % FLINT_BITS; len2 = limbs2 - s2;
   m = FLINT_MIN(s1, s2);
   
   /* this many output limbs will be zero */
   flint_mpn_zero(arrayg, m);

   /* set in1 to shifted array1 */
   if (temp != NULL)
      in1 = temp;
   else
      in1 = flint_malloc(len1*sizeof(mp_limb_t));
   if (b1 == 0)
      flint_mpn_copyi(in1, array1 + s1, len1);
   else
      mpn_rshift(in1, array1 + s1, len1, b1);
   len1 -= (in1[len1 - 1] == 0); 

   /* set in2 to shifted array2 */
   if (temp != NULL)
      in2 = temp + len1;
   else
      in2 = flint_malloc(len2*sizeof(mp_limb_t));
   if (b2 == 0)
      flint_mpn_copyi(in2, array2 + s2, len2);
   else
      mpn_rshift(in2, array2 + s2, len2, b2);
   len2 -= (in2[len2 - 1] == 0); 
   
   
   /* compute gcd of shifted values */
   if (len1 >= len2)
      leng = mpn_gcd(arrayg + m, in1, len1, in2, len2);
   else 
      leng = mpn_gcd(arrayg + m, in2, len2, in1, len1);

   if (mb) /* shift back by mb bits */
   {
      cy = mpn_lshift(arrayg + m, arrayg + m, leng, mb);
      if (cy)
         arrayg[m + leng++] = cy;
   }

   /* clean up */
   if (temp == NULL)
   {
      flint_free(in1);
      flint_free(in2);
   }

   /* return total number of limbs in output */
   return m + leng;
}

mp_size_t flint_mpn_gcd_full(mp_ptr arrayg,
        mp_srcptr array1, mp_size_t limbs1, mp_srcptr array2, mp_size_t limbs2)
{
   return flint_mpn_gcd_full2(arrayg, array1, limbs1, array2, limbs2, NULL);
}

