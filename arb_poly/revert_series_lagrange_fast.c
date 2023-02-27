/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

/* pointer to (x/Q)^i */
#define Ri(ii) (R + (n-1)*((ii)-1))

void
_arb_poly_revert_series_lagrange_fast(arb_ptr Qinv, arb_srcptr Q, slong Qlen, slong n, slong prec)
{
    slong i, j, m;
    arb_ptr R, S, T, tmp;
    arb_t t;

    if (n <= 2)
    {
        if (n >= 1)
            arb_zero(Qinv);
        if (n == 2)
            arb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    m = n_sqrt(n);

    arb_init(t);
    R = _arb_vec_init((n - 1) * m);
    S = _arb_vec_init(n - 1);
    T = _arb_vec_init(n - 1);

    arb_zero(Qinv);
    arb_inv(Qinv + 1, Q + 1, prec);

    _arb_poly_inv_series(Ri(1), Q + 1, FLINT_MIN(Qlen, n) - 1, n - 1, prec);
    for (i = 2; i <= m; i++)
        _arb_poly_mullow(Ri(i), Ri((i + 1) / 2), n - 1, Ri(i / 2), n - 1, n - 1, prec);

    for (i = 2; i < m; i++)
        arb_div_ui(Qinv + i, Ri(i) + i - 1, i, prec);

    _arb_vec_set(S, Ri(m), n - 1);

    for (i = m; i < n; i += m)
    {
        arb_div_ui(Qinv + i, S + i - 1, i, prec);

        for (j = 1; j < m && i + j < n; j++)
        {
            arb_dot(t, NULL, 0, S, 1, Ri(j) + i + j - 1, -1, i + j, prec);
            arb_div_ui(Qinv + i + j, t, i + j, prec);
        }

        if (i + 1 < n)
        {
            _arb_poly_mullow(T, S, n - 1, Ri(m), n - 1, n - 1, prec);
            tmp = S; S = T; T = tmp;
        }
    }

    arb_clear(t);
    _arb_vec_clear(R, (n - 1) * m);
    _arb_vec_clear(S, n - 1);
    _arb_vec_clear(T, n - 1);
}

void
arb_poly_revert_series_lagrange_fast(arb_poly_t Qinv,
                                    const arb_poly_t Q, slong n, slong prec)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !arb_is_zero(Q->coeffs)
                 || arb_contains_zero(Q->coeffs + 1))
    {
        flint_printf("Exception (arb_poly_revert_series_lagrange_fast). Input \n"
               "must have zero constant term and nonzero coefficient of x^1.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        arb_poly_fit_length(Qinv, n);
        _arb_poly_revert_series_lagrange_fast(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_revert_series_lagrange_fast(t->coeffs, Q->coeffs, Qlen, n, prec);
        arb_poly_swap(Qinv, t);
        arb_poly_clear(t);
    }

    _arb_poly_set_length(Qinv, n);
    _arb_poly_normalise(Qinv);
}

