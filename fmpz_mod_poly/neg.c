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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void _fmpz_mod_poly_neg(fmpz *res, const fmpz *poly, len_t len, const fmpz_t p)
{
    len_t i;

    for (i = 0; i < len; i++)
    {
        if (poly[i])
            fmpz_sub(res + i, p, poly + i);
        else
            fmpz_zero(res + i);
    }
}

void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly)
{
    const len_t len = poly->length;

    fmpz_mod_poly_fit_length(res, len);
    _fmpz_mod_poly_set_length(res, len);

    _fmpz_mod_poly_neg(res->coeffs, poly->coeffs, poly->length, &(poly->p));
}

