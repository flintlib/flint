/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2008, 2009 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

// Assumes poly1 and poly2 are not length 0 and len1 >= len2
void _fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, ulong len1, 
					                  const fmpz * poly2, ulong len2)
{
   int neg1, neg2;
   ulong i, limbs1, limbs2, log_len;
   long bits1 = 0, bits2 = 0, bits;
   mp_limb_t * arr1, * arr2, * arr3;
   long sign = 0;
   ulong coeff;
   int sgn;
   
   for (coeff = len1 - 1; coeff >= 0; coeff--)
      if ((neg1 = fmpz_sgn(poly1 + coeff)) != 0) break;
   if (neg1 >= 0) neg1 = 0;

   for (coeff = len2 - 1; coeff >= 0; coeff--)
      if ((neg2 = fmpz_sgn(poly2 + coeff)) != 0) break;
   if (neg2 >= 0) neg2 = 0;

   bits1 = _fmpz_vec_max_bits(poly1, len1);
   if (bits1 < 0L) { sign = 1; bits1 = -bits1; }
 
   if (poly1 != poly2)
   {
	  bits2 = _fmpz_vec_max_bits(poly2, len2);
      if (bits2 < 0L) { sign = 1; bits2 = -bits2; }
   } else
	  bits2 = bits1;
   
   log_len = FLINT_BIT_COUNT(len2);
   bits = bits1 + bits2 + log_len + sign;

   limbs1 = (bits*len1 - 1)/FLINT_BITS + 1;
   limbs2 = (bits*len2 - 1)/FLINT_BITS + 1;
   
   arr1 = (mp_limb_t *) calloc(limbs1, sizeof(mp_limb_t));
   _fmpz_poly_bit_pack(arr1, poly1, len1, bits, neg1);
   
   if (poly1 != poly2)
   {
	  arr2 = (mp_limb_t *) calloc(limbs2, sizeof(mp_limb_t));
      _fmpz_poly_bit_pack(arr2, poly2, len2, bits, neg2);
   }

   arr3 = (mp_limb_t *) malloc((limbs1 + limbs2)*sizeof(mp_limb_t));

   if (poly1 != poly2)
	  mpn_mul(arr3, arr1, limbs1, arr2, limbs2);
   else
	  mpn_mul_n(arr3, arr1, arr1, limbs1);

   if (sign)
	  _fmpz_poly_bit_unpack(res, len1 + len2 - 1, arr3, bits, neg1 ^ neg2);
   else
	  _fmpz_poly_bit_unpack_unsigned(res, len1 + len2 - 1, arr3, bits);

   free(arr1);
   if (poly1 != poly2) free(arr2);
   free(arr3);
}

void fmpz_poly_mul_KS(fmpz_poly_t res, 
                         const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   ulong len_out;
   
   if ((poly1->length == 0) || (poly2->length == 0))  
   {
      fmpz_poly_zero(res);
      return;
   }

   len_out = poly1->length + poly2->length - 1;

   if (res == poly1 || res == poly2)
   {
	   fmpz_poly_t temp;
	   fmpz_poly_init(temp);
	   fmpz_poly_fit_length(temp, len_out);
       if (poly1->length >= poly2->length)
		  _fmpz_poly_mul_KS(temp->coeffs, poly1->coeffs, poly1->length,
		                              poly2->coeffs, poly2->length);
	   else
		  _fmpz_poly_mul_KS(temp->coeffs, poly2->coeffs, poly2->length,
		                              poly1->coeffs, poly1->length);
	   fmpz_poly_swap(res, temp);
	   fmpz_poly_clear(temp);
   } else
   {
	   fmpz_poly_fit_length(res, len_out);
       if (poly1->length >= poly2->length)
		  _fmpz_poly_mul_KS(res->coeffs, poly1->coeffs, poly1->length,
		                              poly2->coeffs, poly2->length);
	   else
		  _fmpz_poly_mul_KS(res->coeffs, poly2->coeffs, poly2->length,
		                              poly1->coeffs, poly1->length);
   }

   _fmpz_poly_set_length(res, len_out);
   _fmpz_poly_normalise(res);
}
