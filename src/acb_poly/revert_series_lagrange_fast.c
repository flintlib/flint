/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

/* pointer to (x/Q)^i */
#define Ri(ii) (R + (n-1)*((ii)-1))

void
_acb_poly_revert_series_lagrange_fast(acb_ptr Qinv, acb_srcptr Q, slong Qlen, slong n, slong prec)
{
    slong i, j, m;
    acb_ptr R, S, T, tmp;
    acb_t t;

    if (n <= 2)
    {
        if (n >= 1)
            acb_zero(Qinv);
        if (n == 2)
            acb_inv(Qinv + 1, Q + 1, prec);
        return;
    }

    m = n_sqrt(n);

    acb_init(t);
    R = _acb_vec_init((n - 1) * m);
    S = _acb_vec_init(n - 1);
    T = _acb_vec_init(n - 1);

    acb_zero(Qinv);
    acb_inv(Qinv + 1, Q + 1, prec);

    _acb_poly_inv_series(Ri(1), Q + 1, FLINT_MIN(Qlen, n) - 1, n - 1, prec);
    for (i = 2; i <= m; i++)
        _acb_poly_mullow(Ri(i), Ri((i + 1) / 2), n - 1, Ri(i / 2), n - 1, n - 1, prec);

    for (i = 2; i < m; i++)
        acb_div_ui(Qinv + i, Ri(i) + i - 1, i, prec);

    _acb_vec_set(S, Ri(m), n - 1);

    for (i = m; i < n; i += m)
    {
        acb_div_ui(Qinv + i, S + i - 1, i, prec);

        for (j = 1; j < m && i + j < n; j++)
        {
            acb_dot(t, NULL, 0, S, 1, Ri(j) + i + j - 1, -1, i + j, prec);
            acb_div_ui(Qinv + i + j, t, i + j, prec);
        }

        if (i + 1 < n)
        {
            _acb_poly_mullow(T, S, n - 1, Ri(m), n - 1, n - 1, prec);
            tmp = S; S = T; T = tmp;
        }
    }

    acb_clear(t);
    _acb_vec_clear(R, (n - 1) * m);
    _acb_vec_clear(S, n - 1);
    _acb_vec_clear(T, n - 1);
}

void
acb_poly_revert_series_lagrange_fast(acb_poly_t Qinv,
                                    const acb_poly_t Q, slong n, slong prec)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !acb_is_zero(Q->coeffs)
                 || acb_contains_zero(Q->coeffs + 1))
    {
        flint_printf("Exception (acb_poly_revert_series_lagrange_fast). Input \n"
               "must have zero constant term and nonzero coefficient of x^1.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        acb_poly_fit_length(Qinv, n);
        _acb_poly_revert_series_lagrange_fast(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_revert_series_lagrange_fast(t->coeffs, Q->coeffs, Qlen, n, prec);
        acb_poly_swap(Qinv, t);
        acb_poly_clear(t);
    }

    _acb_poly_set_length(Qinv, n);
    _acb_poly_normalise(Qinv);
}

