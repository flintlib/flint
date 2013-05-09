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

    Copyright (C) 2011, 2010 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, len_t len)
{
    len_t i;

    for (i = 1; i < len; i++)
        fmpz_mul_ui(rpoly + (i - 1), poly + i, i);
}

void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly)
{
    const len_t len = poly->length;

    if (len < 2)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        fmpz_poly_fit_length(res, len - 1);
        _fmpz_poly_derivative(res->coeffs, poly->coeffs, len);
        _fmpz_poly_set_length(res, len - 1);
    }
}
