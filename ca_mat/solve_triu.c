/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_solve_triu_classical(ca_mat_t X, const ca_mat_t U,
    const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    slong i, j, n, m;
    ca_ptr tmp;
    ca_t s;

    n = U->r;
    m = B->c;

    ca_init(s, ctx);
    tmp = flint_malloc(sizeof(ca_struct) * n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = *ca_mat_entry(X, j, i);

        for (j = n - 1; j >= 0; j--)
        {
            ca_dot(s, ca_mat_entry(B, j, i), 1, U->rows[j] + j + 1, 1, tmp + j + 1, 1, n - j - 1, ctx);

            if (!unit)
                ca_div(tmp + j, s, ca_mat_entry(U, j, j), ctx);
            else
                ca_swap(tmp + j, s, ctx);
        }

        for (j = 0; j < n; j++)
            *ca_mat_entry(X, j, i) = tmp[j];
    }

    flint_free(tmp);
    ca_clear(s, ctx);
}

void
ca_mat_solve_triu_recursive(ca_mat_t X,
        const ca_mat_t U, const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    ca_mat_t UA, UB, UD, XX, XY, BX, BY, T;
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

    ca_mat_window_init(UA, U, 0, 0, r, r, ctx);
    ca_mat_window_init(UB, U, 0, r, r, n, ctx);
    ca_mat_window_init(UD, U, r, r, n, n, ctx);
    ca_mat_window_init(BX, B, 0, 0, r, m, ctx);
    ca_mat_window_init(BY, B, r, 0, n, m, ctx);
    ca_mat_window_init(XX, X, 0, 0, r, m, ctx);
    ca_mat_window_init(XY, X, r, 0, n, m, ctx);

    ca_mat_solve_triu(XY, UD, BY, unit, ctx);

    /* ca_mat_submul(XX, BX, UB, XY); */
    ca_mat_init(T, UB->r, XY->c, ctx);
    ca_mat_mul(T, UB, XY, ctx);
    ca_mat_sub(XX, BX, T, ctx);
    ca_mat_clear(T, ctx);

    ca_mat_solve_triu(XX, UA, XX, unit, ctx);

    ca_mat_window_clear(UA, ctx);
    ca_mat_window_clear(UB, ctx);
    ca_mat_window_clear(UD, ctx);
    ca_mat_window_clear(BX, ctx);
    ca_mat_window_clear(BY, ctx);
    ca_mat_window_clear(XX, ctx);
    ca_mat_window_clear(XY, ctx);
}

void
ca_mat_solve_triu(ca_mat_t X, const ca_mat_t U,
                                    const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    /* todo: tune thresholds */
    if (B->r < 10 || B->c < 10)
        ca_mat_solve_triu_classical(X, U, B, unit, ctx);
    else
        ca_mat_solve_triu_recursive(X, U, B, unit, ctx);
}
