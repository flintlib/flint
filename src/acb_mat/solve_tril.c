/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_solve_tril_classical(acb_mat_t X,
        const acb_mat_t L, const acb_mat_t B, int unit, slong prec)
{
    slong i, j, n, m;
    acb_ptr tmp;
    acb_t s;

    n = L->r;
    m = B->c;

    acb_init(s);
    tmp = flint_malloc(sizeof(acb_struct) * n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = *acb_mat_entry(X, j, i);

        for (j = 0; j < n; j++)
        {
            acb_dot(s, acb_mat_entry(B, j, i), 1, L->rows[j], 1, tmp, 1, j, prec);

            if (!unit)
                acb_div(tmp + j, s, acb_mat_entry(L, j, j), prec);
            else
                acb_swap(tmp + j, s);
        }

        for (j = 0; j < n; j++)
            *acb_mat_entry(X, j, i) = tmp[j];
    }

    flint_free(tmp);
    acb_clear(s);
}

void
acb_mat_solve_tril_recursive(acb_mat_t X,
        const acb_mat_t L, const acb_mat_t B, int unit, slong prec)
{
    acb_mat_t LA, LC, LD, XX, XY, BX, BY, T;
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
    acb_mat_window_init(LA, L, 0, 0, r, r);
    acb_mat_window_init(LC, L, r, 0, n, r);
    acb_mat_window_init(LD, L, r, r, n, n);
    acb_mat_window_init(BX, B, 0, 0, r, m);
    acb_mat_window_init(BY, B, r, 0, n, m);
    acb_mat_window_init(XX, X, 0, 0, r, m);
    acb_mat_window_init(XY, X, r, 0, n, m);

    acb_mat_solve_tril(XX, LA, BX, unit, prec);

    /* acb_mat_submul(XY, BY, LC, XX); */
    acb_mat_init(T, LC->r, BX->c);
    acb_mat_mul(T, LC, XX, prec);
    acb_mat_sub(XY, BY, T, prec);
    acb_mat_clear(T);

    acb_mat_solve_tril(XY, LD, XY, unit, prec);

    acb_mat_window_clear(LA);
    acb_mat_window_clear(LC);
    acb_mat_window_clear(LD);
    acb_mat_window_clear(BX);
    acb_mat_window_clear(BY);
    acb_mat_window_clear(XX);
    acb_mat_window_clear(XY);
}

void
acb_mat_solve_tril(acb_mat_t X, const acb_mat_t L,
                                    const acb_mat_t B, int unit, slong prec)
{
    if (B->r < 40 || B->c < 40)
        acb_mat_solve_tril_classical(X, L, B, unit, prec);
    else
        acb_mat_solve_tril_recursive(X, L, B, unit, prec);
}
