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
ca_mat_solve_tril_classical(ca_mat_t X,
        const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    slong i, j, n, m;
    ca_ptr tmp;
    ca_t s;

    n = L->r;
    m = B->c;

    ca_init(s, ctx);
    tmp = flint_malloc(sizeof(ca_struct) * n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = *ca_mat_entry(X, j, i);

        for (j = 0; j < n; j++)
        {
            ca_dot(s, ca_mat_entry(B, j, i), 1, L->rows[j], 1, tmp, 1, j, ctx);

            if (!unit)
                ca_div(tmp + j, s, ca_mat_entry(L, j, j), ctx);
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
ca_mat_solve_tril_recursive(ca_mat_t X,
        const ca_mat_t L, const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    ca_mat_t LA, LC, LD, XX, XY, BX, BY, T;
    slong r, n, m;

    n = L->r;
    m = B->c;
    r = n / 2;

    if (n == 0 || m == 0)
        return;

    /*
    Denoting inv(M) by M^, we have:

    [A 0]^ [X]  ==  [A^          0 ] [X]  ==  [A^ X]
    [C D]  [Y]  ==  [-D^ C A^    D^] [Y]  ==  [D^ (Y - C A^ X)]
    */
    ca_mat_window_init(LA, L, 0, 0, r, r, ctx);
    ca_mat_window_init(LC, L, r, 0, n, r, ctx);
    ca_mat_window_init(LD, L, r, r, n, n, ctx);
    ca_mat_window_init(BX, B, 0, 0, r, m, ctx);
    ca_mat_window_init(BY, B, r, 0, n, m, ctx);
    ca_mat_window_init(XX, X, 0, 0, r, m, ctx);
    ca_mat_window_init(XY, X, r, 0, n, m, ctx);

    ca_mat_solve_tril(XX, LA, BX, unit, ctx);

    /* ca_mat_submul(XY, BY, LC, XX); */
    ca_mat_init(T, LC->r, BX->c, ctx);
    ca_mat_mul(T, LC, XX, ctx);
    ca_mat_sub(XY, BY, T, ctx);
    ca_mat_clear(T, ctx);

    ca_mat_solve_tril(XY, LD, XY, unit, ctx);

    ca_mat_window_clear(LA, ctx);
    ca_mat_window_clear(LC, ctx);
    ca_mat_window_clear(LD, ctx);
    ca_mat_window_clear(BX, ctx);
    ca_mat_window_clear(BY, ctx);
    ca_mat_window_clear(XX, ctx);
    ca_mat_window_clear(XY, ctx);
}

void
ca_mat_solve_tril(ca_mat_t X, const ca_mat_t L,
                                    const ca_mat_t B, int unit, ca_ctx_t ctx)
{
    /* todo: tune thresholds */
    if (B->r < 10 || B->c < 10)
        ca_mat_solve_tril_classical(X, L, B, unit, ctx);
    else
        ca_mat_solve_tril_recursive(X, L, B, unit, ctx);
}
