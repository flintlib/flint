/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

static void
_apply_permutation(slong * AP, ca_mat_t A, slong * P, slong n, slong offset)
{
    if (n != 0)
    {
        ca_struct ** Atmp;
        slong *APtmp;
        slong i;

        Atmp = flint_malloc(sizeof(ca_mat_struct *) * n);
        APtmp = flint_malloc(sizeof(slong) * n);

        for (i = 0; i < n; i++)
            Atmp[i] = A->rows[P[i] + offset];
        for (i = 0; i < n; i++)
            A->rows[i + offset] = Atmp[i];

        for (i = 0; i < n; i++)
            APtmp[i] = AP[P[i] + offset];
        for (i = 0; i < n; i++)
            AP[i + offset] = APtmp[i];

        flint_free(Atmp);
        flint_free(APtmp);
    }
}

int
ca_mat_lu_recursive(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
{

    slong i, j, m, n, r1, r2, n1;
    ca_mat_t A0, A1, A00, A01, A10, A11;
    slong *P1;
    int success;

    r2 = 0;
    m = A->r;
    n = A->c;

    if (m < 4 || n < 4)
        return ca_mat_lu_classical(rank, P, LU, A, rank_check, ctx);

    if (LU != A)
        ca_mat_set(LU, A, ctx);

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(slong) * m);
    ca_mat_window_init(A0, LU, 0, 0, m, n1, ctx);
    ca_mat_window_init(A1, LU, 0, n1, m, n, ctx);

    success = ca_mat_lu_recursive(&r1, P1, A0, A0, rank_check, ctx);

    if (!success || (rank_check && (r1 != n1)))
    {
        success = 0;
        goto cleanup1;
    }

    if (r1 != 0)
    {
        _apply_permutation(P, LU, P1, m, 0);
    }

    ca_mat_window_init(A00, LU, 0, 0, r1, r1, ctx);
    ca_mat_window_init(A10, LU, r1, 0, m, r1, ctx);
    ca_mat_window_init(A01, LU, 0, n1, r1, n, ctx);
    ca_mat_window_init(A11, LU, r1, n1, m, n, ctx);

    if (r1 != 0)
    {
        ca_mat_t T;
        ca_mat_init(T, ca_mat_nrows(A10), ca_mat_ncols(A01), ctx);
        ca_mat_solve_tril(A01, A00, A01, 1, ctx);
        ca_mat_mul(T, A10, A01, ctx);
        ca_mat_sub(A11, A11, T, ctx);
        ca_mat_clear(T, ctx);
    }

    success = ca_mat_lu_recursive(&r2, P1, A11, A11, rank_check, ctx);

    if (!success || (rank_check && (r1 + r2 < FLINT_MIN(m, n))))
    {
        r1 = r2 = 0;
        goto cleanup2;
    }

    _apply_permutation(P, LU, P1, m - r1, r1);

    /* Compress L */
    if (r1 != n1)
    {
        for (i = 0; i < m - r1; i++)
        {
            ca_struct * row = LU->rows[r1 + i];

            for (j = 0; j < FLINT_MIN(i, r2) ; j++)
            {
                ca_set(row + r1 + j, row + n1 + j, ctx);
                ca_zero(row + n1 + j, ctx);
            }
        }
    }

cleanup2:
    ca_mat_window_clear(A00, ctx);
    ca_mat_window_clear(A10, ctx);
    ca_mat_window_clear(A01, ctx);
    ca_mat_window_clear(A11, ctx);

cleanup1:
    flint_free(P1);
    ca_mat_window_clear(A0, ctx);
    ca_mat_window_clear(A1, ctx);

    *rank = r1 + r2;

    return success;
}
