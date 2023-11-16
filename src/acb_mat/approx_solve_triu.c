/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
acb_approx_mul(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)),
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARB_RND);
}

/* note: the tmp variable t should have zero radius */
static void
acb_approx_div(acb_t z, const acb_t x, const acb_t y, acb_t t, slong prec)
{
    arf_set(arb_midref(acb_realref(t)), arb_midref(acb_realref(y)));
    arf_set(arb_midref(acb_imagref(t)), arb_midref(acb_imagref(y)));

    acb_inv(t, t, prec);

    mag_zero(arb_radref(acb_realref(t)));
    mag_zero(arb_radref(acb_imagref(t)));

    acb_approx_mul(z, x, t, prec);
}

void
acb_mat_approx_solve_triu_classical(acb_mat_t X, const acb_mat_t U,
    const acb_mat_t B, int unit, slong prec)
{
    slong i, j, n, m;
    acb_ptr tmp;
    acb_t s, t;

    n = U->r;
    m = B->c;

    acb_init(s);
    acb_init(t);
    tmp = flint_malloc(sizeof(acb_struct) * n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = *acb_mat_entry(X, j, i);

        for (j = n - 1; j >= 0; j--)
        {
            acb_approx_dot(s, acb_mat_entry(B, j, i), 1, U->rows[j] + j + 1, 1, tmp + j + 1, 1, n - j - 1, prec);

            if (!unit)
                acb_approx_div(tmp + j, s, arb_mat_entry(U, j, j), t, prec);
            else
                acb_swap(tmp + j, s);
        }

        for (j = 0; j < n; j++)
            *acb_mat_entry(X, j, i) = tmp[j];
    }

    flint_free(tmp);
    acb_clear(s);
    acb_clear(t);
}

void
acb_mat_approx_solve_triu_recursive(acb_mat_t X,
        const acb_mat_t U, const acb_mat_t B, int unit, slong prec)
{
    acb_mat_t UA, UB, UD, XX, XY, BX, BY, T;
    slong r, n, m;

    n = U->r;
    m = B->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    /*
    Denoting inv(M) by M^, we have:
    [A B]^ [X]  ==  [A^ (X - B D^ Y)]
    [0 D]  [Y]  ==  [    D^ Y       ]
    */

    acb_mat_window_init(UA, U, 0, 0, r, r);
    acb_mat_window_init(UB, U, 0, r, r, n);
    acb_mat_window_init(UD, U, r, r, n, n);
    acb_mat_window_init(BX, B, 0, 0, r, m);
    acb_mat_window_init(BY, B, r, 0, n, m);
    acb_mat_window_init(XX, X, 0, 0, r, m);
    acb_mat_window_init(XY, X, r, 0, n, m);

    acb_mat_approx_solve_triu(XY, UD, BY, unit, prec);

    acb_mat_init(T, UB->r, XY->c);
    acb_mat_approx_mul(T, UB, XY, prec);
    acb_mat_sub(XX, BX, T, prec);
    acb_mat_get_mid(XX, XX);
    acb_mat_clear(T);

    acb_mat_approx_solve_triu(XX, UA, XX, unit, prec);

    acb_mat_window_clear(UA);
    acb_mat_window_clear(UB);
    acb_mat_window_clear(UD);
    acb_mat_window_clear(BX);
    acb_mat_window_clear(BY);
    acb_mat_window_clear(XX);
    acb_mat_window_clear(XY);
}

void
acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U,
                                    const acb_mat_t B, int unit, slong prec)
{
    if (B->r < 40 || B->c < 40)
        acb_mat_approx_solve_triu_classical(X, U, B, unit, prec);
    else
        acb_mat_approx_solve_triu_recursive(X, U, B, unit, prec);
}
