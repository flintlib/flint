/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include <string.h>
#include "templates.h"


static void
_apply_permutation(slong * AP, TEMPLATE(T, mat_t) A, const slong * P,
    slong num_rows, slong row_offset, slong num_cols, slong col_offset)
{
    if (num_rows != 0)
    {
        TEMPLATE(T, struct) * Atmp;
        slong * APtmp;
        slong i;

        /* todo: reduce memory allocation */
        Atmp = flint_malloc(sizeof(TEMPLATE(T, struct)) * num_rows * num_cols);
        /* todo: avoid temporary allocation when AP != P */
        APtmp = flint_malloc(sizeof(slong) * num_rows);

        for (i = 0; i < num_rows; i++)
            memcpy(Atmp + i * num_cols, TEMPLATE(T, mat_entry) (A, P[i] + row_offset, col_offset), num_cols * sizeof(TEMPLATE(T, struct)));
        for (i = 0; i < num_rows; i++)
            memcpy(TEMPLATE(T, mat_entry) (A, i + row_offset, col_offset), Atmp + i * num_cols, num_cols * sizeof(TEMPLATE(T, struct)));

        for (i = 0; i < num_rows; i++) APtmp[i] = AP[P[i] + row_offset];
        for (i = 0; i < num_rows; i++) AP[i + row_offset] = APtmp[i];

        flint_free(Atmp);
        flint_free(APtmp);
    }
}


slong
TEMPLATE(T, mat_lu_recursive) (slong * P,
                               TEMPLATE(T, mat_t) A,
                               int rank_check, const TEMPLATE(T, ctx_t) ctx)
{

    slong i, j, m, n, r1, r2, n1;
    TEMPLATE(T, mat_t) A0, A1, A00, A01, A10, A11;
    slong *P1;

    m = A->r;
    n = A->c;

    if (m < TEMPLATE(CAP_T, MAT_LU_RECURSIVE_CUTOFF)
        || n < TEMPLATE(CAP_T, MAT_LU_RECURSIVE_CUTOFF))
    {
        r1 = TEMPLATE(T, mat_lu_classical) (P, A, rank_check, ctx);
        return r1;
    }

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(slong) * m);
    TEMPLATE(T, mat_window_init) (A0, A, 0, 0, m, n1, ctx);
    TEMPLATE(T, mat_window_init) (A1, A, 0, n1, m, n, ctx);

    r1 = TEMPLATE(T, mat_lu) (P1, A0, rank_check, ctx);

    if (rank_check && (r1 != n1))
    {
        flint_free(P1);
        TEMPLATE(T, mat_window_clear) (A0, ctx);
        TEMPLATE(T, mat_window_clear) (A1, ctx);
        return 0;
    }

    if (r1 != 0)
    {
        _apply_permutation(P, A, P1, m, 0, n - n1, n1);
    }

    TEMPLATE(T, mat_window_init) (A00, A, 0, 0, r1, r1, ctx);
    TEMPLATE(T, mat_window_init) (A10, A, r1, 0, m, r1, ctx);
    TEMPLATE(T, mat_window_init) (A01, A, 0, n1, r1, n, ctx);
    TEMPLATE(T, mat_window_init) (A11, A, r1, n1, m, n, ctx);

    if (r1 != 0)
    {
        TEMPLATE(T, mat_solve_tril) (A01, A00, A01, 1, ctx);
        TEMPLATE(T, mat_submul) (A11, A11, A10, A01, ctx);
    }

    r2 = TEMPLATE(T, mat_lu) (P1, A11, rank_check, ctx);

    if (rank_check && (r1 + r2 < FLINT_MIN(m, n)))
    {
        r1 = r2 = 0;
    }
    else
    {
        _apply_permutation(P, A, P1, m - r1, r1, n1, 0);

        /* Compress L */
        if (r1 != n1)
        {
            for (i = 0; i < m - r1; i++)
            {
                TEMPLATE(T, struct) * row = TEMPLATE(T, mat_entry) (A, r1 + i, 0);
                for (j = 0; j < FLINT_MIN(i, r2); j++)
                {
                    TEMPLATE(T, set) (row + r1 + j, row + n1 + j, ctx);
                    TEMPLATE(T, zero) (row + n1 + j, ctx);
                }
            }
        }
    }

    flint_free(P1);
    TEMPLATE(T, mat_window_clear) (A00, ctx);
    TEMPLATE(T, mat_window_clear) (A01, ctx);
    TEMPLATE(T, mat_window_clear) (A10, ctx);
    TEMPLATE(T, mat_window_clear) (A11, ctx);
    TEMPLATE(T, mat_window_clear) (A0, ctx);
    TEMPLATE(T, mat_window_clear) (A1, ctx);

    return r1 + r2;
}


#endif
