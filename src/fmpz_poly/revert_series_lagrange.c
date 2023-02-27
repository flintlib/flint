/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_revert_series_lagrange(fmpz * Qinv,
    const fmpz * Q, slong Qlen, slong n)
{
    slong i;
    fmpz *R, *S, *T, *tmp;

    if (n <= 2)
    {
        _fmpz_vec_set(Qinv, Q, n);
        return;
    }

    R = _fmpz_vec_init(n - 1);
    S = _fmpz_vec_init(n - 1);
    T = _fmpz_vec_init(n - 1);

    fmpz_zero(Qinv);
    fmpz_set(Qinv + 1, Q + 1);

    _fmpz_poly_inv_series(R, Q + 1, FLINT_MIN(Qlen, n) - 1, n - 1);
    _fmpz_vec_set(S, R, n - 1);

    for (i = 2; i < n; i++)
    {
        _fmpz_poly_mullow(T, S, n - 1, R, n - 1, n - 1);
        fmpz_divexact_ui(Qinv + i, T + i - 1, i);
        tmp = S; S = T; T = tmp;
    }

    _fmpz_vec_clear(R, n - 1);
    _fmpz_vec_clear(S, n - 1);
    _fmpz_vec_clear(T, n - 1);
}

void
fmpz_poly_revert_series_lagrange(fmpz_poly_t Qinv,
                                        const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !fmpz_is_zero(Q->coeffs) || !fmpz_is_pm1(Q->coeffs + 1))
    {
        flint_printf("Exception (fmpz_poly_revert_series_lagrange). Input must have \n"
               "zero constant term and +1 or -1 as coefficient of x^1\n.");
        flint_abort();
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_revert_series_lagrange(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_revert_series_lagrange(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }
    
    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}
