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
#include "fmpz_poly.h"

void _fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, long len)
{
    fmpz * high = res + len;
    fmpz_t x;
    fmpz_init(x);
    _fmpz_poly_content(x, poly, len);
    if (fmpz_sgn(poly + (len - 1)) < 0)
        fmpz_neg(x, x);
    if (*x != 1L)
        while (res != high)
            fmpz_divexact(res++, poly++, x);
    else if (res != poly)
        while (res != high)
            fmpz_set(res++, poly++);
    fmpz_clear(x);
}

void fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly)
{
    long len = poly->length;
    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    
    fmpz_poly_fit_length(res, len);
    _fmpz_poly_set_length(res, len);
    
    _fmpz_poly_primitive_part(res->coeffs, poly->coeffs, len);
}

