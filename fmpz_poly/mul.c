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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, 
				   ulong len1, const fmpz * poly2, ulong len2)
{
   if (len2 == 1)
   {
      _fmpz_poly_mul_classical(res, poly1, len1, poly2, len2);
	  return;
   }

   ulong len_out = len1 + len2 - 1;
   ulong limbs1 = _fmpz_vec_max_limbs(poly1, len1);
   ulong limbs2 = _fmpz_vec_max_limbs(poly2, len2);
   ulong max_limbs = FLINT_MAX(limbs1, limbs2);

   if ((max_limbs > 4) && (len_out < 16))
      _fmpz_poly_mul_karatsuba(res, poly1, len1, poly2, len2);
   else
      _fmpz_poly_mul_KS(res, poly1, len1, poly2, len2);

}

void fmpz_poly_mul(fmpz_poly_t res, 
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
		  _fmpz_poly_mul(temp->coeffs, poly1->coeffs, poly1->length,
		                              poly2->coeffs, poly2->length);
	   else
		  _fmpz_poly_mul(temp->coeffs, poly2->coeffs, poly2->length,
		                              poly1->coeffs, poly1->length);
	   fmpz_poly_swap(res, temp);
	   fmpz_poly_clear(temp);
   } else
   {
	   fmpz_poly_fit_length(res, len_out);
       if (poly1->length >= poly2->length)
		  _fmpz_poly_mul(res->coeffs, poly1->coeffs, poly1->length,
		                              poly2->coeffs, poly2->length);
	   else
		  _fmpz_poly_mul(res->coeffs, poly2->coeffs, poly2->length,
		                              poly1->coeffs, poly1->length);
   }

   _fmpz_poly_set_length(res, len_out);
   _fmpz_poly_normalise(res);
}
