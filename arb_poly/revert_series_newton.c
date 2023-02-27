/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

#define CUTOFF 5

void
_arb_poly_revert_series_newton(arb_ptr Qinv, arb_srcptr Q, slong Qlen, slong n, slong prec)
{
    slong i, k, a[FLINT_BITS];
    arb_ptr T, U, V;

    if (n <= 2)
    {
        if (n >= 1)
            arb_zero(Qinv);
        if (n == 2)
            arb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    T = _arb_vec_init(n);
    U = _arb_vec_init(n);
    V = _arb_vec_init(n);

    k = n;
    for (i = 1; (WORD(1) << i) < k; i++);
    a[i = 0] = k;
    while (k >= CUTOFF)
        a[++i] = (k = (k + 1) / 2);

    _arb_poly_revert_series_lagrange(Qinv, Q, Qlen, k, prec);
    _arb_vec_zero(Qinv + k, n - k);

    for (i--; i >= 0; i--)
    {
        k = a[i];
        _arb_poly_compose_series(T, Q, FLINT_MIN(Qlen, k), Qinv, k, k, prec);
        _arb_poly_derivative(U, T, k, prec); arb_zero(U + k - 1);
        arb_zero(T + 1);
        _arb_poly_div_series(V, T, k, U, k, k, prec);
        _arb_poly_derivative(T, Qinv, k, prec);
        _arb_poly_mullow(U, V, k, T, k, k, prec);
        _arb_vec_sub(Qinv, Qinv, U, k, prec);
    }

    _arb_vec_clear(T, n);
    _arb_vec_clear(U, n);
    _arb_vec_clear(V, n);
}

void
arb_poly_revert_series_newton(arb_poly_t Qinv,
                                    const arb_poly_t Q, slong n, slong prec)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !arb_is_zero(Q->coeffs)
                 || arb_contains_zero(Q->coeffs + 1))
    {
        flint_printf("Exception (arb_poly_revert_series_newton). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        arb_poly_fit_length(Qinv, n);
        _arb_poly_revert_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_revert_series_newton(t->coeffs, Q->coeffs, Qlen, n, prec);
        arb_poly_swap(Qinv, t);
        arb_poly_clear(t);
    }

    _arb_poly_set_length(Qinv, n);
    _arb_poly_normalise(Qinv);
}

