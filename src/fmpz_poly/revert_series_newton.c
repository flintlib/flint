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

#define FLINT_REVERSE_NEWTON_CUTOFF 10

void
_fmpz_poly_revert_series_newton(fmpz * Qinv,
    const fmpz * Q, slong Qlen, slong n)
{
    fmpz *T, *U, *V;
    slong alloc = 3 * n;

    if (n <= 2)
    {
        _fmpz_vec_set(Qinv, Q, n);
        return;
    }

    T = _fmpz_vec_init(alloc);
    U = T + n;
    V = U + n;

    FLINT_NEWTON_INIT(FLINT_REVERSE_NEWTON_CUTOFF, n)

    FLINT_NEWTON_BASECASE(k)
    _fmpz_poly_revert_series_lagrange(Qinv, Q, Qlen, k);
    _fmpz_vec_zero(Qinv + k, n - k);
    FLINT_NEWTON_END_BASECASE

    FLINT_NEWTON_LOOP(FLINT_UNUSED(k0), k)
    _fmpz_poly_compose_series(T, Q, FLINT_MIN(Qlen, k), Qinv, k, k);
    _fmpz_poly_derivative(U, T, k); fmpz_zero(U + k - 1);
    fmpz_zero(T + 1);
    _fmpz_poly_div_series(V, T, k, U, k, k);
    _fmpz_poly_derivative(T, Qinv, k);
    _fmpz_poly_mullow(U, V, k, T, k, k);
    _fmpz_vec_sub(Qinv, Qinv, U, k);
    FLINT_NEWTON_END_LOOP

    FLINT_NEWTON_END

    _fmpz_vec_clear(T, alloc);
}

void
fmpz_poly_revert_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !fmpz_is_zero(Q->coeffs) || !fmpz_is_pm1(Q->coeffs + 1))
    {
        flint_printf("Exception (fmpz_poly_revert_series_newton). Input must have \n"
               "zero constant term and +1 or -1 as coefficient of x^1\n.");
        flint_abort();
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_revert_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_revert_series_newton(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }
    
    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}
