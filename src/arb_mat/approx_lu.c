/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

static void
_apply_permutation(slong * AP, arb_mat_t A, slong * P,
    slong n, slong offset)
{
    if (n != 0)
    {
        arb_ptr * Atmp;
        slong * APtmp;
        slong i;

        Atmp = flint_malloc(sizeof(arb_ptr) * n);
        APtmp = flint_malloc(sizeof(slong) * n);

        for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i] + offset];
        for (i = 0; i < n; i++) A->rows[i + offset] = Atmp[i];

        for (i = 0; i < n; i++) APtmp[i] = AP[P[i] + offset];
        for (i = 0; i < n; i++) AP[i + offset] = APtmp[i];

        flint_free(Atmp);
        flint_free(APtmp);
    }
}

void
_arb_vec_approx_scalar_addmul(arb_ptr res, arb_srcptr vec,
    slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arf_addmul(arb_midref(res + i),
                   arb_midref(vec + i), arb_midref(c), prec, ARF_RND_DOWN);
}

int
arb_mat_approx_lu_classical(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    arf_t d;
    arb_t e;
    arb_ptr * a;
    slong i, j, m, n, r, row, col;
    int result;

    if (arb_mat_is_empty(A))
        return 1;

    m = arb_mat_nrows(A);
    n = arb_mat_ncols(A);

    arb_mat_get_mid(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    arf_init(d);
    arb_init(e);

    result = 1;

    while (row < m && col < n)
    {
        r = arb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            arb_mat_swap_rows(LU, P, row, r);

        arf_ui_div(d, 1, arb_midref(a[row] + col), prec, ARB_RND);

        for (j = row + 1; j < m; j++)
        {
            arf_mul(arb_midref(e), arb_midref(a[j] + col), d, prec, ARB_RND);
            arb_neg(e, e);
            _arb_vec_approx_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            arf_zero(arb_midref(a[j] + col));
            arb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    arf_clear(d);
    arb_clear(e);

    return result;
}

int
arb_mat_approx_lu_recursive(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    slong i, m, n, r1, r2, n1;
    arb_mat_t A0, A1, A00, A01, A10, A11;
    slong * P1;

    m = A->r;
    n = A->c;

    if (m <= 1 || n <= 1)
    {
        return arb_mat_approx_lu_classical(P, LU, A, prec);
    }

    arb_mat_get_mid(LU, A);

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(slong) * m);
    arb_mat_window_init(A0, LU, 0, 0, m, n1);
    arb_mat_window_init(A1, LU, 0, n1, m, n);

    r1 = arb_mat_approx_lu(P1, A0, A0, prec);

    if (!r1)
    {
        flint_free(P1);
        arb_mat_window_clear(A0);
        arb_mat_window_clear(A1);
        return 0;
    }

    /* r1 = rank of A0 */
    r1 = FLINT_MIN(m, n1);

    _apply_permutation(P, LU, P1, m, 0);

    arb_mat_window_init(A00, LU, 0, 0, r1, r1);
    arb_mat_window_init(A10, LU, r1, 0, m, r1);
    arb_mat_window_init(A01, LU, 0, n1, r1, n);
    arb_mat_window_init(A11, LU, r1, n1, m, n);

    arb_mat_approx_solve_tril(A01, A00, A01, 1, prec);

    {
        /* arb_mat_approx_submul(A11, A11, A10, A01, prec); */
        arb_mat_t T;
        arb_mat_init(T, A10->r, A01->c);
        arb_mat_approx_mul(T, A10, A01, prec);
        arb_mat_sub(A11, A11, T, prec);
        arb_mat_get_mid(A11, A11);
        arb_mat_clear(T);
    }

    r2 = arb_mat_approx_lu(P1, A11, A11, prec);

    if (!r2)
        r1 = r2 = 0;
    else
        _apply_permutation(P, LU, P1, m - r1, r1);

    flint_free(P1);
    arb_mat_window_clear(A00);
    arb_mat_window_clear(A01);
    arb_mat_window_clear(A10);
    arb_mat_window_clear(A11);
    arb_mat_window_clear(A0);
    arb_mat_window_clear(A1);

    return r1 && r2;
}

int
arb_mat_approx_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    if (arb_mat_nrows(A) < 8 || arb_mat_ncols(A) < 8)
        return arb_mat_approx_lu_classical(P, LU, A, prec);
    else
        return arb_mat_approx_lu_recursive(P, LU, A, prec);
}
