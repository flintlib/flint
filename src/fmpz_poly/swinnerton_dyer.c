/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "arb_poly.h"

void
_fmpz_poly_swinnerton_dyer(fmpz * T, ulong n)
{
    if (n == 0)
    {
        fmpz_zero(T);
        fmpz_one(T + 1);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_swinnerton_dyer_ui(t, n, 0);
        if (!_arb_vec_get_unique_fmpz_vec(T, t->coeffs, t->length))
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);
        arb_poly_clear(t);
    }
}

void
fmpz_poly_swinnerton_dyer(fmpz_poly_t poly, ulong n)
{
    slong N = (WORD(1) << n);
    fmpz_poly_fit_length(poly, N + 1);
    _fmpz_poly_swinnerton_dyer(poly->coeffs, n);
    _fmpz_poly_set_length(poly, N + 1);
}

