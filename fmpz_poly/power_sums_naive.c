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

    Copyright (C) 2016 Vincent Delecroix

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_power_sums_naive(fmpz * res, const fmpz * poly, slong len, slong n)
{
    slong i, k;

    fmpz_set_ui(res, len - 1);
    for (k = 1; k < FLINT_MIN(n, len); k++)
    {
        fmpz_mul_ui(res + k, poly + len - 1 - k, k);
        for (i = 1; i < k; i++)
            fmpz_addmul(res + k, poly + len - 1 - k + i, res + i);
        fmpz_neg(res + k, res + k);
    }
    for (k = len; k < n; k++)
    {
        fmpz_zero(res + k);
        for (i = k - len + 1; i < k; i++)
            fmpz_addmul(res + k, poly + len - 1 - k + i, res + i);
        fmpz_neg(res + k, res + k);
    }
}

void
fmpz_poly_power_sums_naive(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    if (poly->length == 0)
    {
        flint_printf
            ("Exception (fmpz_poly_power_sums_naive). Zero polynomial.\n");
        flint_abort();
    }
    else if (n <= 0 || poly->length == 1)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        if (poly == res)
        {
            fmpz_poly_t t;
            fmpz_poly_init(t);
            fmpz_poly_fit_length(t, n);
            _fmpz_poly_power_sums_naive(t->coeffs, poly->coeffs,
                                        poly->length, n);
            fmpz_poly_swap(res, t);
            fmpz_poly_clear(t);
        }
        else
        {
            fmpz_poly_fit_length(res, n);
            _fmpz_poly_power_sums_naive(res->coeffs, poly->coeffs,
                                        poly->length, n);
        }
        _fmpz_poly_set_length(res, n);
        _fmpz_poly_normalise(res);
    }
}
