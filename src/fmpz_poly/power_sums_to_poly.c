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
_fmpz_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, slong len)
{
    slong k;
    slong d = fmpz_get_ui(poly);

    fmpz_one(res + d);
    for (k = 1; k < FLINT_MIN(d + 1, len); k++)
    {
        _fmpz_vec_dot_general(res + d - k, poly + k, 0, res + d - k + 1, poly + 1, 0, k - 1);
        fmpz_divexact_si(res + d - k, res + d - k, -k);
    }
    for (k = len; k <= d; k++)
    {
        _fmpz_vec_dot_general(res + d - k, NULL, 0, res + d - k + 1, poly + 1, 0, len - 1);
        fmpz_divexact_si(res + d - k, res + d - k, -k);
    }
}

void
fmpz_poly_power_sums_to_poly(fmpz_poly_t res, const fmpz_poly_t Q)
{
    if (Q->length == 0)
    {
        fmpz_poly_fit_length(res, 1);
        fmpz_one(res->coeffs);
        _fmpz_poly_set_length(res, 1);
    }
    else
    {
        slong d;
        d = fmpz_get_ui(Q->coeffs);
        if (Q == res)
        {
            fmpz_poly_t t;
            fmpz_poly_init(t);
            fmpz_poly_fit_length(t, d + 1);
            _fmpz_poly_power_sums_to_poly(t->coeffs, Q->coeffs, Q->length);
            fmpz_poly_swap(res, t);
            fmpz_poly_clear(t);
        }
        else
        {
            fmpz_poly_fit_length(res, d + 1);
            _fmpz_poly_power_sums_to_poly(res->coeffs, Q->coeffs, Q->length);
        }
        _fmpz_poly_set_length(res, d + 1);
        _fmpz_poly_normalise(res);
    }
}
