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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_reverse(fmpz * res, const fmpz * poly, len_t len, len_t n)
{
    if (res == poly)
    {
        len_t i;

        for (i = 0; i < n / 2; i++)
        {
            fmpz t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            fmpz_zero(res + i);
    }
    else
    {
        len_t i;

        for (i = 0; i < n - len; i++)
            fmpz_zero(res + i);

        for (i = 0; i < len; i++)
            fmpz_set(res + (n - len) + i, poly + (len - 1) - i);
    }
}

void
fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, len_t n)
{
    len_t len = FLINT_MIN(n, poly->length);
    if (len == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, n);

    _fmpz_poly_reverse(res->coeffs, poly->coeffs, len, n);

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
