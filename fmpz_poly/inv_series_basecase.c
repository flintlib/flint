/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_inv_series_basecase(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
{
    Qlen = FLINT_MIN(Qlen, n);

    fmpz_set(Qinv, Q);

    if (Qlen == 1)
    {
        _fmpz_vec_zero(Qinv + 1, n - 1);
    }
    else
    {
        slong i, j;

        for (i = 1; i < n; i++)
        {
            fmpz_mul(Qinv + i, Q + 1, Qinv + i - 1);

            for (j = 2; j < FLINT_MIN(i + 1, Qlen); j++)
                fmpz_addmul(Qinv + i, Q + j, Qinv + i - j);

            if (fmpz_is_one(Qinv))
                fmpz_neg(Qinv + i, Qinv + i);
        }
    }
}

void
fmpz_poly_inv_series_basecase(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_printf("Exception (fmpz_poly_inv_series_basecase). Division by zero.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_inv_series_basecase(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_inv_series_basecase(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }

    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}

