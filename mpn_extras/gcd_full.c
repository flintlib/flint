/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

mp_size_t flint_mpn_gcd_full(mp_ptr arrayg, 
    mp_ptr array1, mp_size_t limbs1, mp_ptr array2, mp_size_t limbs2)
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
   if (b1 == 0)
      in1 = array1 + s1;
   else
   {
      in1 = flint_malloc(len1*sizeof(mp_limb_t));
      mpn_rshift(in1, array1 + s1, len1, b1);
      len1 -= (in1[len1 - 1] == 0); 
   }

   /* set in2 to shifted array2 */
   if (b2 == 0)
      in2 = array2 + s2;
   else
   {
      in2 = flint_malloc(len2*sizeof(mp_limb_t));
      mpn_rshift(in2, array2 + s2, len2, b2);
      len2 -= (in2[len2 - 1] == 0); 
   }
   
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
   if (b1) flint_free(in1);
   if (b2) flint_free(in2);

   /* return total number of limbs in output */
   return m + leng;
}
