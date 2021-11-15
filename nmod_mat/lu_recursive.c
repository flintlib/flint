/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


static void
_apply_permutation(slong * AP, nmod_mat_t A, slong * P,
    slong n, slong offset)
{
    if (n != 0)
    {
        mp_ptr * Atmp;
        slong * APtmp;
        slong i;

        Atmp = flint_malloc(sizeof(mp_ptr) * n);
        APtmp = flint_malloc(sizeof(slong) * n);

        for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i] + offset];
        for (i = 0; i < n; i++) A->rows[i + offset] = Atmp[i];

        for (i = 0; i < n; i++) APtmp[i] = AP[P[i] + offset];
        for (i = 0; i < n; i++) AP[i + offset] = APtmp[i];

        flint_free(Atmp);
        flint_free(APtmp);
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
        _apply_permutation(P, A, P1, m, 0);
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
        _apply_permutation(P, A, P1, m - r1, r1);

        /* Compress L */
        if (r1 != n1)
        {
            for (i = 0; i < m - r1; i++)
            {
                mp_ptr row = A->rows[r1 + i];
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
