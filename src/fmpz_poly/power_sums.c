/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_power_sums(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    slong len = poly->length;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_power_sums): Zero polynomial.\n");
    }
    else if (n <= 0 || len == 1)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        size_t i = 0;
        while (fmpz_is_zero(poly->coeffs + i))
            i++;
        if (poly == res)
        {
            fmpz_poly_t t;
            fmpz_poly_init2(t, n);
            _fmpz_poly_power_sums_naive(t->coeffs, poly->coeffs + i,
                                           len - i, n);
            fmpz_poly_swap(res, t);
            fmpz_poly_clear(t);
        }
        else
        {
            fmpz_poly_fit_length(res, n);
            _fmpz_poly_power_sums_naive(res->coeffs, poly->coeffs + i,
                                           len - i, n);
        }
        _fmpz_poly_set_length(res, n);
        if (i)
            fmpz_set_si(res->coeffs, len - 1);
        _fmpz_poly_normalise(res);
    }
}
