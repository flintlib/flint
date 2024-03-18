/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_power_sums_naive(fmpz * res, const fmpz * poly, slong len, slong n)
{
    slong k;

    fmpz_set_ui(res, len - 1);

    for (k = 1; k < FLINT_MIN(n, len); k++)
    {
        fmpz_mul_si(res + k, poly + len - 1 - k, -k);
        _fmpz_vec_dot_general(res + k, res + k, 1, poly + len - 1 - k + 1, res + 1, 0, k - 1);
    }

    for (k = len; k < n; k++)
        _fmpz_vec_dot_general(res + k, NULL, 1, poly, res + k - len + 1, 0, len - 1);
}

void
fmpz_poly_power_sums_naive(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_power_sums_naive): Zero polynomial.\n");
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
