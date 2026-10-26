/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "nmod_mat/impl.h"

static slong
_ple_sub(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots, slong cutoff)
{
    if (cutoff <= 0)
        return _nmod_mat_ple(P, Q, A, max_nonpivots);
    else if (FLINT_MIN(A->r, A->c) >= cutoff)
        return _nmod_mat_ple_recursive_cutoff(P, Q, A, max_nonpivots, cutoff);
    else
        return _nmod_mat_ple_classical_delayed(P, Q, A, max_nonpivots);
}

/*
    The column-split recursion of nmod_mat_lu_recursive. The left half is
    decomposed, the right half gets the same row order, its first r1 rows
    are solved with the unit lower triangular L of the left half and its
    other rows updated with the rest of that L, and the bottom right block
    is decomposed. Its row order is applied to the L of the left half, and
    its L moves to the columns r1, r1 + 1, ... The rows without pivot of
    both decompositions are in increasing order, hence so are those of the
    result: the output is the canonical decomposition of nmod_mat_ple.
*/
slong
_nmod_mat_ple_recursive_cutoff(slong * P, slong * Q, nmod_mat_t A,
                               slong max_nonpivots, slong cutoff)
{
    const slong m = A->r, n = A->c;
    slong n1, r1, r2, i, j, k, max_nonpivots2;
    slong * P2, * Q2, * T;
    nmod_mat_t A0, A1, A00, A01, A10, A11;
    TMP_INIT;

    if (m <= 1 || n <= 1)
        return _nmod_mat_ple_classical(P, Q, A, max_nonpivots);

    n1 = n / 2;

    /* left half */
    nmod_mat_window_init(A0, A, 0, 0, m, n1);
    r1 = _ple_sub(P, Q, A0, max_nonpivots, cutoff);
    nmod_mat_window_clear(A0);

    if (r1 < 0)
        return -1;

    /* right half: rows in the same order */
    nmod_mat_window_init(A1, A, 0, n1, m, n);
    _nmod_mat_ple_permute_rows(A1, P);
    nmod_mat_window_clear(A1);

    nmod_mat_window_init(A00, A, 0, 0, r1, r1);
    nmod_mat_window_init(A10, A, r1, 0, m, r1);
    nmod_mat_window_init(A01, A, 0, n1, r1, n);
    nmod_mat_window_init(A11, A, r1, n1, m, n);

    if (r1 != 0)
    {
        nmod_mat_solve_tril(A01, A00, A01, 1);
        nmod_mat_submul(A11, A11, A10, A01);
    }

    /* remaining bound: the left half had n1 - r1 columns without pivot,
       all met while rows remained (r1 < m); if r1 = m, no row remains */
    if (max_nonpivots < 0 || r1 == m)
        max_nonpivots2 = -1;
    else
        max_nonpivots2 = max_nonpivots - (n1 - r1);

    TMP_START;
    P2 = TMP_ALLOC((m - r1 + 1) * sizeof(slong));
    Q2 = TMP_ALLOC((n - n1 + 1) * sizeof(slong));
    T = TMP_ALLOC((FLINT_MAX(m, n) + 1) * sizeof(slong));

    r2 = _ple_sub(P2, Q2, A11, max_nonpivots2, cutoff);

    if (r2 >= 0)
    {
        /* rows r1, r1 + 1, ... of the left half (L) in the new order */
        _nmod_mat_ple_permute_rows(A10, P2);
        for (i = 0; i < m - r1; i++)
            T[i] = P[r1 + P2[i]];
        for (i = 0; i < m - r1; i++)
            P[r1 + i] = T[i];

        /* move the L of the bottom right block to the columns r1, ... */
        if (r1 != n1)
        {
            for (i = 1; i < m - r1; i++)
            {
                ulong * row = nmod_mat_entry_ptr(A, r1 + i, 0);

                for (j = 0; j < FLINT_MIN(i, r2); j++)
                {
                    row[r1 + j] = row[n1 + j];
                    row[n1 + j] = 0;
                }
            }
        }

        /* pivot columns, then the others, in increasing order */
        for (k = 0; k < n1 - r1; k++)
            T[k] = Q[r1 + k];
        for (k = 0; k < r2; k++)
            Q[r1 + k] = n1 + Q2[k];
        for (k = 0; k < n1 - r1; k++)
            Q[r1 + r2 + k] = T[k];
        for (k = r2; k < n - n1; k++)
            Q[n1 + k] = n1 + Q2[k];
    }

    TMP_END;

    nmod_mat_window_clear(A00);
    nmod_mat_window_clear(A01);
    nmod_mat_window_clear(A10);
    nmod_mat_window_clear(A11);

    return (r2 < 0) ? -1 : r1 + r2;
}

slong
_nmod_mat_ple_recursive(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    return _nmod_mat_ple_recursive_cutoff(P, Q, A, max_nonpivots, 0);
}

slong
nmod_mat_ple_recursive(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
{
    return _nmod_mat_ple_wrap(P, pivots_nonpivots, A, flags,
                              _nmod_mat_ple_recursive);
}
