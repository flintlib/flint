/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

static void
_apply_permutation(slong * AP, gr_mat_t A, slong * P, slong n, slong offset)
{
    if (n != 0)
    {
        gr_ptr * Atmp;
        slong *APtmp;
        slong i;

        /* todo: stack alloc */
        Atmp = flint_malloc(sizeof(gr_ptr) * n);
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
gr_mat_lu_recursive(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong i, j, m, n, r1, r2, n1;
    gr_mat_t A0, A1, A00, A01, A10, A11;
    slong *P1;
    int status = GR_SUCCESS;

    m = A->r;
    n = A->c;

    if (m < 4 || n < 4)
        return gr_mat_lu_classical(rank, P, LU, A, rank_check, ctx);

    if (LU != A)
        status |= gr_mat_set(LU, A, ctx);

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    r1 = r2 = 0;

    P1 = flint_malloc(sizeof(slong) * m);
    gr_mat_window_init(A0, LU, 0, 0, m, n1, ctx);
    gr_mat_window_init(A1, LU, 0, n1, m, n, ctx);

    status |= gr_mat_lu_recursive(&r1, P1, A0, A0, rank_check, ctx);

    if (status != GR_SUCCESS)
        goto cleanup1;

    /* We proved that the matrix is rank-deficient, accomplishing the goal. */
    if (rank_check && r1 != n1)
    {
        r1 = r2 = 0;
        goto cleanup1;
    }

    if (r1 != 0)
    {
        _apply_permutation(P, LU, P1, m, 0);
    }

    gr_mat_window_init(A00, LU, 0, 0, r1, r1, ctx);
    gr_mat_window_init(A10, LU, r1, 0, m, r1, ctx);
    gr_mat_window_init(A01, LU, 0, n1, r1, n, ctx);
    gr_mat_window_init(A11, LU, r1, n1, m, n, ctx);

    if (r1 != 0)
    {
        gr_mat_t T;
        gr_mat_init(T, A10->r, A01->c, ctx);
        status |= gr_mat_nonsingular_solve_tril(A01, A00, A01, 1, ctx);
        status |= gr_mat_mul(T, A10, A01, ctx);
        status |= gr_mat_sub(A11, A11, T, ctx);
        gr_mat_clear(T, ctx);
    }

    status |= gr_mat_lu_recursive(&r2, P1, A11, A11, rank_check, ctx);

    if (status != GR_SUCCESS)
        goto cleanup2;

    /* We proved that the matrix is rank-deficient, accomplishing the goal. */
    if (rank_check && (r1 + r2 < FLINT_MIN(m, n)))
    {
        r1 = r2 = 0;
        goto cleanup2;
    }

    _apply_permutation(P, LU, P1, m - r1, r1);

    /* Compress L */
    if (r1 != n1)
    {
        slong sz = ctx->sizeof_elem;

        for (i = 0; i < m - r1; i++)
        {
            gr_ptr row = LU->rows[r1 + i];

            for (j = 0; j < FLINT_MIN(i, r2); j++)
            {
                status |= gr_set(GR_ENTRY(row, r1 + j, sz), GR_ENTRY(row, n1 + j, sz), ctx);
                status |= gr_zero(GR_ENTRY(row, n1 + j, sz), ctx);
            }
        }
    }

cleanup2:
    gr_mat_window_clear(A00, ctx);
    gr_mat_window_clear(A10, ctx);
    gr_mat_window_clear(A01, ctx);
    gr_mat_window_clear(A11, ctx);

cleanup1:
    flint_free(P1);
    gr_mat_window_clear(A0, ctx);
    gr_mat_window_clear(A1, ctx);

    *rank = r1 + r2;

    return status;
}
