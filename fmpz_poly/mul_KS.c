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

// Assumes poly1 and poly2 are not length 0
void _fmpz_poly_mul_KS(fmpz_poly_t res, 
                         const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   int neg1, neg2;
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;
   ulong i, limbs1, limbs2, log_len;
   mp_bitcnt_t bits1 = 0, bits2 = 0, bits;
   mp_limb_t * arr1, * arr2, * arr3;

   neg1 = fmpz_sgn(poly1->coeffs + len1 - 1);
   if (neg1 > 0) neg1 = 0;

   for (i = 0; i < len1; i++)
      bits1 = FLINT_MAX(bits1, fmpz_bits(poly1->coeffs + i));
   
   neg2 = fmpz_sgn(poly2->coeffs + len2 - 1);
   if (neg2 > 0) neg2 = 0;

   for (i = 0; i < len2; i++)
      bits2 = FLINT_MAX(bits2, fmpz_bits(poly2->coeffs + i));
   
   log_len = FLINT_BIT_COUNT(len2);
   bits = bits1 + bits2 + log_len + 1;

   limbs1 = (bits*len1 - 1)/FLINT_BITS + 1;
   limbs2 = (bits*len2 - 1)/FLINT_BITS + 1;
   arr1 = (mp_limb_t *) calloc(limbs1, sizeof(mp_limb_t));
   arr2 = (mp_limb_t *) calloc(limbs2, sizeof(mp_limb_t));
   
   arr3 = (mp_limb_t *) malloc((limbs1 + limbs2)*sizeof(mp_limb_t));

   fmpz_poly_bit_pack(arr1, poly1, bits, neg1);
   fmpz_poly_bit_pack(arr2, poly2, bits, neg2);

   mpn_mul(arr3, arr1, limbs1, arr2, limbs2);

   fmpz_poly_bit_unpack(res, len1 + len2 - 1, arr3, bits, neg1 ^ neg2);

   free(arr1);
   free(arr2);
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
		  _fmpz_poly_mul_KS(temp, poly1, poly2);
	   else
		  _fmpz_poly_mul_KS(temp, poly2, poly1);
	   fmpz_poly_swap(res, temp);
	   fmpz_poly_clear(temp);
   } else
   {
	   fmpz_poly_fit_length(res, len_out);
       if (poly1->length >= poly2->length)
		  _fmpz_poly_mul_KS(res, poly1, poly2);
	   else
		  _fmpz_poly_mul_KS(res, poly2, poly1);
   }

   _fmpz_poly_set_length(res, len_out);
   _fmpz_poly_normalise(res);
}
