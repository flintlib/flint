/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

    Loosely based on the recursive PLS implementation in M4RI,
    Copyright (C) 2008 Clement Pernet.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


static void
_apply_permutation(len_t * AP, nmod_mat_t A, len_t * P,
    len_t n, len_t offset)
{
    if (n != 0)
    {
        mp_ptr * Atmp;
        len_t * APtmp;
        len_t i;

        Atmp = flint_malloc(sizeof(mp_ptr) * n);
        APtmp = flint_malloc(sizeof(len_t) * n);

        for (i = 0; i < n; i++) Atmp[i] = A->rows[P[i] + offset];
        for (i = 0; i < n; i++) A->rows[i + offset] = Atmp[i];

        for (i = 0; i < n; i++) APtmp[i] = AP[P[i] + offset];
        for (i = 0; i < n; i++) AP[i + offset] = APtmp[i];

        flint_free(Atmp);
        flint_free(APtmp);
    }
}


len_t 
nmod_mat_lu_recursive(len_t * P, nmod_mat_t A, int rank_check)
{
    len_t i, j, m, n, r1, r2, n1;
    nmod_mat_t A0, A1, A00, A01, A10, A11;
    len_t * P1;

    m = A->r;
    n = A->c;

    if (m < NMOD_MAT_LU_RECURSIVE_CUTOFF || n < NMOD_MAT_LU_RECURSIVE_CUTOFF)
    {
        r1 = nmod_mat_lu_classical(P, A, rank_check);
        return r1;
    }

    n1 = n / 2;

    for (i = 0; i < m; i++)
        P[i] = i;

    P1 = flint_malloc(sizeof(len_t) * m);
    nmod_mat_window_init(A0, A, 0, 0, m, n1);
    nmod_mat_window_init(A1, A, 0, n1, m, n);

    r1 = nmod_mat_lu(P1, A0, rank_check);

    if (rank_check && (r1 != n1))
    {
        flint_free(P1);
        nmod_mat_window_clear(A0);
        nmod_mat_window_clear(A1);
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
    nmod_mat_window_clear(A1);

    return r1 + r2;
}
