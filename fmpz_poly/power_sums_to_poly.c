/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"


void
_fmpz_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, slong len)
{
    slong i, k;
    slong d = fmpz_get_ui(poly);

    fmpz_one(res + d);
    for (k = 1; k < FLINT_MIN(d + 1, len); k++)
    {
        fmpz_set(res + d - k, poly + k);
        for (i = 1; i < k; i++)
            fmpz_addmul(res + d - k, res + d - k + i, poly + i);
        fmpz_divexact_si(res + d - k, res + d - k, k);
        fmpz_neg(res + d - k, res + d - k);
    }
    for (k = len; k <= d; k++)
    {
        fmpz_zero(res + d - k);
        for (i = 1; i < len; i++)
            fmpz_addmul(res + d - k, res + d - k + i, poly + i);
        fmpz_divexact_si(res + d - k, res + d - k, k);
        fmpz_neg(res + d - k, res + d - k);
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
