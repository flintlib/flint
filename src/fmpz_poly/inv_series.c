/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
{
    if (Qlen < 64 || n < 64)
        _fmpz_poly_inv_series_basecase(Qinv, Q, Qlen, n);
    else
        _fmpz_poly_inv_series_newton(Qinv, Q, Qlen, n);
}

void
fmpz_poly_inv_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_inv_series). Division by zero.\n");
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_inv_series(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_inv_series(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }

    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}

