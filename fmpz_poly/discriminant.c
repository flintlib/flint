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

    Copyright (C) 2014 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_discriminant(fmpz_t res, const fmpz * poly, slong len)
{
   fmpz * der = _fmpz_vec_init(len - 1);

   _fmpz_poly_derivative(der, poly, len);
   _fmpz_poly_resultant(res, poly, len, der, len - 1);
   
   if ((len & 3) == 0 || (len & 3) == 3) /* degree is not 0, 1 mod 4 */
      fmpz_neg(res, res);

   if (!fmpz_is_one(poly + len - 1))
      fmpz_divexact(res, res, poly + len - 1);

   _fmpz_vec_clear(der, len - 1);
}

void fmpz_poly_discriminant(fmpz_t res, const fmpz_poly_t poly)
{
   slong len = poly->length;
   
   if (len <= 1)
     fmpz_zero(res);
   else
      _fmpz_poly_discriminant(res, poly->coeffs, len);  
}
