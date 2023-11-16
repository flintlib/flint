/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_solve_tril_classical(arb_mat_t X,
        const arb_mat_t L, const arb_mat_t B, int unit, slong prec)
{
    slong i, j, n, m;
    arb_ptr tmp;
    arb_t s;

    n = L->r;
    m = B->c;

    arb_init(s);
    tmp = flint_malloc(sizeof(arb_struct) * n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = *arb_mat_entry(X, j, i);

        for (j = 0; j < n; j++)
        {
            arb_dot(s, arb_mat_entry(B, j, i), 1, L->rows[j], 1, tmp, 1, j, prec);

            if (!unit)
                arb_div(tmp + j, s, arb_mat_entry(L, j, j), prec);
            else
                arb_swap(tmp + j, s);
        }

        for (j = 0; j < n; j++)
            *arb_mat_entry(X, j, i) = tmp[j];
    }

    flint_free(tmp);
    arb_clear(s);
}

void
arb_mat_solve_tril_recursive(arb_mat_t X,
        const arb_mat_t L, const arb_mat_t B, int unit, slong prec)
{
    arb_mat_t LA, LC, LD, XX, XY, BX, BY, T;
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
    arb_mat_window_init(LA, L, 0, 0, r, r);
    arb_mat_window_init(LC, L, r, 0, n, r);
    arb_mat_window_init(LD, L, r, r, n, n);
    arb_mat_window_init(BX, B, 0, 0, r, m);
    arb_mat_window_init(BY, B, r, 0, n, m);
    arb_mat_window_init(XX, X, 0, 0, r, m);
    arb_mat_window_init(XY, X, r, 0, n, m);

    arb_mat_solve_tril(XX, LA, BX, unit, prec);

    /* arb_mat_submul(XY, BY, LC, XX); */
    arb_mat_init(T, LC->r, BX->c);
    arb_mat_mul(T, LC, XX, prec);
    arb_mat_sub(XY, BY, T, prec);
    arb_mat_clear(T);

    arb_mat_solve_tril(XY, LD, XY, unit, prec);

    arb_mat_window_clear(LA);
    arb_mat_window_clear(LC);
    arb_mat_window_clear(LD);
    arb_mat_window_clear(BX);
    arb_mat_window_clear(BY);
    arb_mat_window_clear(XX);
    arb_mat_window_clear(XY);
}

void
arb_mat_solve_tril(arb_mat_t X, const arb_mat_t L,
                                    const arb_mat_t B, int unit, slong prec)
{
    if (B->r < 40 || B->c < 40)
        arb_mat_solve_tril_classical(X, L, B, unit, prec);
    else
        arb_mat_solve_tril_recursive(X, L, B, unit, prec);
}

