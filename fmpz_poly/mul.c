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

void fmpz_poly_mul(fmpz_poly_t res, 
             const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;
   ulong len_out = len1 + len2 - 1;
   ulong limbs1 = fmpz_poly_max_limbs(poly1);
   ulong limbs2 = fmpz_poly_max_limbs(poly2);
   ulong max_limbs = FLINT_MAX(limbs1, limbs2);

   if ((len1 <= 1) || (len2 <= 1))
   {
	   fmpz_poly_mul_classical(res, poly1, poly2);
	   return;
   }
  
   if ((max_limbs > 4) && (len_out < 16))
      fmpz_poly_mul_karatsuba(res, poly1, poly2);
   else
      fmpz_poly_mul_KS(res, poly1, poly2);
}
