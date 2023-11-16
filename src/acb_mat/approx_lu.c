/*
    Copyright (C) 2018 arbguest

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
_apply_permutation(slong * AP, acb_mat_t A, slong * P,
    slong n, slong offset)
{
    if (n != 0)
    {
        acb_ptr * Atmp;
        slong * APtmp;
        slong i;

        Atmp = flint_malloc(sizeof(acb_ptr) * n);
        APtmp = flint_malloc(sizeof(slong) * n);

        for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i] + offset];
        for (i = 0; i < n; i++) A->rows[i + offset] = Atmp[i];

        for (i = 0; i < n; i++) APtmp[i] = AP[P[i] + offset];
        for (i = 0; i < n; i++) AP[i + offset] = APtmp[i];

        flint_free(Atmp);
        flint_free(APtmp);
    }
}

static void
_acb_approx_mul(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)),
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARB_RND);
}

static void
_acb_approx_inv(acb_t z, const acb_t x, slong prec)
{
    arf_set(arb_midref(acb_realref(z)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(z)), arb_midref(acb_imagref(x)));

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));

    acb_inv(z, z, prec);

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));
}

static void
_acb_vec_approx_scalar_addmul(acb_ptr res, acb_srcptr vec,
    slong len, const acb_t c, slong prec)
{
    acb_t t;
    slong i;
    acb_init(t);

    for (i = 0; i < len; i++)
    {
        _acb_approx_mul(t, vec + i, c, prec);

        arf_add(arb_midref(acb_realref(res + i)),
            arb_midref(acb_realref(res + i)),
            arb_midref(acb_realref(t)), prec, ARB_RND);
        arf_add(arb_midref(acb_imagref(res + i)),
            arb_midref(acb_imagref(res + i)),
            arb_midref(acb_imagref(t)), prec, ARB_RND);
    }

    acb_clear(t);
}

int
acb_mat_approx_lu_classical(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    acb_t d, e;
    acb_ptr * a;
    slong i, j, m, n, r, row, col;
    int result;

    if (acb_mat_is_empty(A))
        return 1;

    m = acb_mat_nrows(A);
    n = acb_mat_ncols(A);

    acb_mat_get_mid(LU, A);

    a = LU->rows;

    row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    acb_init(d);
    acb_init(e);

    result = 1;

    while (row < m && col < n)
    {
        r = acb_mat_find_pivot_partial(LU, row, m, col);

        if (r == -1)
        {
            result = 0;
            break;
        }
        else if (r != row)
            acb_mat_swap_rows(LU, P, row, r);

        _acb_approx_inv(d, a[row] + col, prec);

        for (j = row + 1; j < m; j++)
        {
            _acb_approx_mul(e, a[j] + col, d, prec);
            acb_neg(e, e);
            _acb_vec_approx_scalar_addmul(a[j] + col,
                a[row] + col, n - col, e, prec);
            acb_zero(a[j] + col);
            acb_neg(a[j] + row, e);
        }

        row++;
        col++;
    }

    acb_clear(d);
    acb_clear(e);

    return result;
}

int
acb_mat_approx_lu_recursive(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    slong i, m, n, r1, r2, n1;
    acb_mat_t A0, A1, A00, A01, A10, A11;
    slong * P1;

    m = A->r;
    n = A->c;

    if (m <= 1 || n <= 1)
    {
        return acb_mat_approx_lu_classical(P, LU, A, prec);
    }

    acb_mat_get_mid(LU, A);

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(slong) * m);
    acb_mat_window_init(A0, LU, 0, 0, m, n1);
    acb_mat_window_init(A1, LU, 0, n1, m, n);

    r1 = acb_mat_approx_lu(P1, A0, A0, prec);

    if (!r1)
    {
        flint_free(P1);
        acb_mat_window_clear(A0);
        acb_mat_window_clear(A1);
        return 0;
    }

    /* r1 = rank of A0 */
    r1 = FLINT_MIN(m, n1);

    _apply_permutation(P, LU, P1, m, 0);

    acb_mat_window_init(A00, LU, 0, 0, r1, r1);
    acb_mat_window_init(A10, LU, r1, 0, m, r1);
    acb_mat_window_init(A01, LU, 0, n1, r1, n);
    acb_mat_window_init(A11, LU, r1, n1, m, n);

    acb_mat_approx_solve_tril(A01, A00, A01, 1, prec);

    {
        /* acb_mat_submul(A11, A11, A10, A01, prec); */
        acb_mat_t T;
        acb_mat_init(T, A10->r, A01->c);
        acb_mat_approx_mul(T, A10, A01, prec);
        acb_mat_sub(A11, A11, T, prec);
        acb_mat_get_mid(A11, A11);
        acb_mat_clear(T);
    }

    r2 = acb_mat_approx_lu(P1, A11, A11, prec);

    if (!r2)
        r1 = r2 = 0;
    else
        _apply_permutation(P, LU, P1, m - r1, r1);

    flint_free(P1);
    acb_mat_window_clear(A00);
    acb_mat_window_clear(A01);
    acb_mat_window_clear(A10);
    acb_mat_window_clear(A11);
    acb_mat_window_clear(A0);
    acb_mat_window_clear(A1);

    return r1 && r2;
}

int
acb_mat_approx_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    if (acb_mat_nrows(A) < 8 || acb_mat_ncols(A) < 8)
        return acb_mat_approx_lu_classical(P, LU, A, prec);
    else
        return acb_mat_approx_lu_recursive(P, LU, A, prec);
}
