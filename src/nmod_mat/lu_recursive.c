/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

/* Permute rows SLICE columns at a time to improve memory locality
   and reduce temporary allocation. */
#define SLICE 32

static void
_apply_permutation(slong * AP, nmod_mat_t A, const slong * P,
    slong num_rows, slong row_offset, slong num_cols, slong col_offset)
{
    if (num_rows != 0)
    {
        ulong * Atmp;
        slong * APtmp;
        slong i, c, l;
        TMP_INIT;

        TMP_START;

        if (num_cols <= SLICE)
        {
            Atmp = TMP_ALLOC(sizeof(ulong) * num_rows * num_cols);

            for (i = 0; i < num_rows; i++)
                _nmod_vec_set(Atmp + i * num_cols, nmod_mat_entry_ptr(A, P[i] + row_offset, col_offset), num_cols);
            for (i = 0; i < num_rows; i++)
                _nmod_vec_set(nmod_mat_entry_ptr(A, i + row_offset, col_offset), Atmp + i * num_cols, num_cols);
        }
        else
        {
            Atmp = TMP_ALLOC(sizeof(ulong) * num_rows * SLICE);

            for (c = 0; c < num_cols; c += SLICE)
            {
                l = FLINT_MIN(SLICE, num_cols - c);
                for (i = 0; i < num_rows; i++)
                    _nmod_vec_set(Atmp + i * SLICE, nmod_mat_entry_ptr(A, P[i] + row_offset, col_offset + c), l);
                for (i = 0; i < num_rows; i++)
                    _nmod_vec_set(nmod_mat_entry_ptr(A, i + row_offset, col_offset + c), Atmp + i * SLICE, l);
            }
        }

        APtmp = TMP_ALLOC(sizeof(slong) * num_rows);

        for (i = 0; i < num_rows; i++) APtmp[i] = AP[P[i] + row_offset];
        for (i = 0; i < num_rows; i++) AP[i + row_offset] = APtmp[i];

        TMP_END;
    }
}



slong
nmod_mat_lu_recursive(slong * P, nmod_mat_t A, int rank_check)
{
    slong i, j, m, n, r1, r2, n1;
    nmod_mat_t A0, A00, A01, A10, A11;
    slong * P1;

    m = A->r;
    n = A->c;

    /* main cutoffs are in nmod_mat_lu */
    if (m <= 1 || n <= 1)
    {
        r1 = nmod_mat_lu_classical(P, A, rank_check);
        return r1;
    }

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(slong) * m);
    nmod_mat_window_init(A0, A, 0, 0, m, n1);

    r1 = nmod_mat_lu(P1, A0, rank_check);

    if (rank_check && (r1 != n1))
    {
        flint_free(P1);
        nmod_mat_window_clear(A0);
        return 0;
    }

    if (r1 != 0)
    {
        _apply_permutation(P, A, P1, m, 0, n - n1, n1);
    }

    nmod_mat_window_init(A00, A, 0, 0, r1, r1);
    nmod_mat_window_init(A10, A, r1, 0, m, r1);
    nmod_mat_window_init(A01, A, 0, n1, r1, n);
    nmod_mat_window_init(A11, A, r1, n1, m, n);

    if (r1 != 0)
    {
        nmod_mat_solve_tril(A01, A00, A01, 1);
        nmod_mat_submul(A11, A11, A10, A01);
    }

    r2 = nmod_mat_lu(P1, A11, rank_check);

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
                nn_ptr row = nmod_mat_entry_ptr(A, r1 + i, 0);

                for (j = 0; j < FLINT_MIN(i, r2); j++)
                {
                    row[r1 + j] = row[n1 + j];
                    row[n1 + j] = 0;
                }
            }
        }
    }

    flint_free(P1);
    nmod_mat_window_clear(A00);
    nmod_mat_window_clear(A01);
    nmod_mat_window_clear(A10);
    nmod_mat_window_clear(A11);
    nmod_mat_window_clear(A0);

    return r1 + r2;
}
