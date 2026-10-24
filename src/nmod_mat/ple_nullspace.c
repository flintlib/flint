/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

/*
    Kernels from a PLE decomposition A = Pi L E of an m x n matrix of rank
    r, as computed by nmod_mat_ple.

    Right kernel: A x = 0 if and only if E x = 0. With E1 (r x r, upper
    triangular, invertible) and E2 (r x (n - r)) the pivot and other
    columns of E, the kernel has the basis with identity at the columns
    without pivot and C = -E1^-1 E2 at the pivot columns.

    Left kernel: x A = 0 if and only if y L = 0 with y_i = x_{P[i]}. With
    L1 (r x r, unit lower triangular) and L2 ((m - r) x r) the first r and
    the other rows of L, the kernel has the basis with identity at the
    rows P[r], ..., P[m-1] and D = -L2 L1^-1 at the rows P[0], ...,
    P[r-1].
*/

/* C = -E1^-1 E2, r x (n - r) */
static void
_ple_right_block(nmod_mat_t C, const nmod_mat_t LU, const slong * Q, slong r)
{
    slong n = LU->c, i, j, k;
    nmod_mat_t E1;

    /* E2: the entries of E in the columns without pivot; those before the
       pivot of the row are zero (and LU may hold L there) */
    for (i = 0; i < r; i++)
        for (k = 0; k < n - r; k++)
            nmod_mat_entry(C, i, k) = (Q[r + k] > Q[i]) ?
                                        nmod_mat_entry(LU, i, Q[r + k]) : 0;

    /* E1, read on and above the diagonal only */
    for (j = 0; j < r && Q[j] == j; j++)
        ;

    if (j == r)
    {
        nmod_mat_window_init(E1, LU, 0, 0, r, r);
        nmod_mat_solve_triu(C, E1, C, 0);
        nmod_mat_window_clear(E1);
    }
    else
    {
        nmod_mat_init(E1, r, r, LU->mod.n);
        for (i = 0; i < r; i++)
            for (j = i; j < r; j++)
                nmod_mat_entry(E1, i, j) = nmod_mat_entry(LU, i, Q[j]);
        nmod_mat_solve_triu(C, E1, C, 0);
        nmod_mat_clear(E1);
    }

    nmod_mat_neg(C, C);
}

/* D = -L2 L1^-1, (m - r) x r */
static void
_ple_left_block(nmod_mat_t D, const nmod_mat_t LU, slong r)
{
    slong m = LU->r;
    nmod_mat_t L1, L2;

    nmod_mat_window_init(L1, LU, 0, 0, r, r);
    nmod_mat_window_init(L2, LU, r, 0, m, r);
    nmod_mat_solve_tril_right(D, L1, L2, 1);
    nmod_mat_neg(D, D);
    nmod_mat_window_clear(L1);
    nmod_mat_window_clear(L2);
}

void
nmod_mat_ple_nullspace_compact(nmod_mat_t X, const nmod_mat_t LU,
                               const slong * pivots_nonpivots, slong rank)
{
    if (rank == 0 || rank == LU->c)
        return;

    _ple_right_block(X, LU, pivots_nonpivots, rank);
}

void
nmod_mat_ple_nullspace(nmod_mat_t X, const nmod_mat_t LU,
                       const slong * pivots_nonpivots, slong rank)
{
    const slong * Q = pivots_nonpivots;
    slong n = LU->c, d = n - rank, j, k;
    nmod_mat_t C;

    if (d == 0)
        return;

    /* the rows of the columns without pivot: identity */
    for (k = 0; k < d; k++)
    {
        _nmod_vec_zero(nmod_mat_entry_ptr(X, Q[rank + k], 0), d);
        nmod_mat_entry(X, Q[rank + k], k) = 1;
    }

    if (rank == 0)
        return;

    nmod_mat_init(C, rank, d, LU->mod.n);
    _ple_right_block(C, LU, Q, rank);
    for (j = 0; j < rank; j++)
        _nmod_vec_set(nmod_mat_entry_ptr(X, Q[j], 0), nmod_mat_entry_ptr(C, j, 0), d);
    nmod_mat_clear(C);
}

void
nmod_mat_ple_left_nullspace_compact(nmod_mat_t X, slong * rows,
                    const nmod_mat_t LU, const slong * P, slong rank)
{
    slong m = LU->r, d = m - rank, i, j, k, t;
    slong * pos;
    nmod_mat_t D;

    /* rows: the pivot rows in increasing order (row rank profile), then
       the others, P[rank..m) */
    pos = flint_malloc((m + 1) * sizeof(slong));
    for (i = 0; i < m; i++)
        pos[i] = -1;
    for (j = 0; j < rank; j++)
        pos[P[j]] = j;
    for (i = t = 0; i < m; i++)
        if (pos[i] >= 0)
            rows[t++] = i;
    for (k = 0; k < d; k++)
        rows[rank + k] = P[rank + k];

    if (rank != 0 && d != 0)
    {
        nmod_mat_init(D, d, rank, LU->mod.n);
        _ple_left_block(D, LU, rank);

        /* column t of X is the column of D for the row rows[t] */
        for (k = 0; k < d; k++)
            for (t = 0; t < rank; t++)
                nmod_mat_entry(X, k, t) = nmod_mat_entry(D, k, pos[rows[t]]);

        nmod_mat_clear(D);
    }

    flint_free(pos);
}

void
nmod_mat_ple_left_nullspace(nmod_mat_t X, const nmod_mat_t LU,
                            const slong * P, slong rank)
{
    slong m = LU->r, d = m - rank, j, k;
    nmod_mat_t D;

    nmod_mat_zero(X);

    if (d == 0)
        return;

    for (k = 0; k < d; k++)
        nmod_mat_entry(X, k, P[rank + k]) = 1;

    if (rank == 0)
        return;

    nmod_mat_init(D, d, rank, LU->mod.n);
    _ple_left_block(D, LU, rank);
    for (k = 0; k < d; k++)
        for (j = 0; j < rank; j++)
            nmod_mat_entry(X, k, P[j]) = nmod_mat_entry(D, k, j);
    nmod_mat_clear(D);
}
