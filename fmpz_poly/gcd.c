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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_gcd(fmpz * res, const fmpz * poly1, long len1, 
                    const fmpz * poly2, long len2)
{
    _fmpz_poly_gcd_subresultant(res, poly1, len1, poly2, len2);
}

void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, 
                   const fmpz_poly_t poly2)
{
    fmpz_poly_gcd_subresultant(res, poly1, poly2);
}
