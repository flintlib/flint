/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"

/*
    Column-wise Crout elimination. The rows that are not pivot rows yet
    are not updated as a whole: when column col is processed, its entries
    in these rows are updated by dot products of their multipliers (L)
    with the column of U, and when a pivot is found, its row is updated in
    all the columns to its right by a vector-matrix product of its
    multipliers with the rows of U. Both products delay their reductions,
    with SIMD code for all moduli once they are long enough (nmod_vec dot
    products, nmod_mat_nmod_vec_mul), and no buffer of unreduced entries
    is needed.

    As in ple_classical_delayed.c, the pivot row k is stored in row k of A
    (the vector-matrix products need the rows of U to form a matrix), a new
    pivot row being swapped with the row at the position of the rank;
    the other rows are stored in any order, orig[i] being the index in the
    input of the row stored in row i, and the pivot is the one with a
    nonzero entry of smallest index in the input (see
    _nmod_mat_ple_rows_finish in impl.h).

    The nonzero multipliers of the row stored in row i are in the columns
    lo[i], ..., hi[i] - 1 (none if lo[i] >= hi[i]). The products of the
    row stop at column hi[i] - 1, and start at column lo[i] when it is at
    least PLE_CROUT_MIN_SKIP: the rows that become zero after a few
    pivots, and the rows with many leading zero multipliers (block or
    staircase structures), cost little. The classical and delayed
    algorithms skip every zero multiplier, which also covers zeros
    scattered among nonzero ones: on matrices with many of these, they do
    less work.
*/

#define PLE_CROUT_MIN_SKIP 16

slong
_nmod_mat_ple_crout(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    const slong m = A->r, n = A->c, stride = A->stride;
    ulong * const aa = A->entries;
    const nmod_t mod = A->mod;
    slong rank = 0, nonpivots = 0, col, i, k, pr;
    slong * orig, * lo, * hi;
    nn_ptr u, tmp;
    dot_params_t params;
    TMP_INIT;

#define ROW(i) (aa + (i) * stride)

    TMP_START;
    orig = TMP_ALLOC(3 * (m + 1) * sizeof(slong));
    lo = orig + m + 1;
    hi = lo + m + 1;
    u = TMP_ALLOC((FLINT_MIN(m, n) + 1) * sizeof(ulong));
    tmp = TMP_ALLOC((n + 1) * sizeof(ulong));
    params = _nmod_vec_dot_params(FLINT_MIN(m, n), mod);

    for (i = 0; i < m; i++)
    {
        orig[i] = i;
        lo[i] = hi[i] = 0;
    }

    for (col = 0; col < n && rank < m; col++)
    {
        ulong * prow;
        ulong inv;

        /* the column in the rows that are not pivot rows, after the
           eliminations of the previous pivots, and the pivot: the row with
           a nonzero entry of smallest index in the input */
        pr = -1;
        if (rank > 0)
        {
            slong best = m;

            for (k = 0; k < rank; k++)
                u[k] = ROW(k)[col];

            for (i = rank; i < m; i++)
            {
                ulong * row = ROW(i);
                ulong x = row[col];
                slong o = orig[i], l = lo[i], h = hi[i];
                int better;

                /* the start of the product is moved only for long runs
                   of zeros: with offset pointers in all cases, this loop
                   was 1.2-1.35x slower on dense matrices, where l = 0
                   (Cascade Lake) */
                if (l >= PLE_CROUT_MIN_SKIP)
                {
                    if (l < h)
                        x = nmod_sub(x, _nmod_vec_dot(row + l, u + l, h - l,
                                                      mod, params), mod);
                }
                else if (h > 0)
                {
                    x = nmod_sub(x, _nmod_vec_dot(row, u, h, mod, params), mod);
                }

                better = (x != 0) & (o < best);

                row[col] = x;
                best = better ? o : best;
                pr = better ? i : pr;
            }
        }
        else
        {
            for (i = 0; i < m && ROW(i)[col] == 0; i++)
                ;
            if (i < m)
                pr = i;
        }

        /* no pivot: the column is zero in the rows that are not pivot rows */
        if (pr < 0)
        {
            nonpivots++;
            if (max_nonpivots >= 0 && nonpivots > max_nonpivots)
            {
                TMP_END;
                return -1;
            }
            continue;
        }

        if (pr != rank)
        {
            nmod_mat_swap_rows(A, NULL, rank, pr);
            FLINT_SWAP(slong, orig[rank], orig[pr]);
            FLINT_SWAP(slong, lo[rank], lo[pr]);
            FLINT_SWAP(slong, hi[rank], hi[pr]);
        }
        P[rank] = orig[rank];
        prow = ROW(rank);

        /* the pivot row to the right of the pivot: minus its multipliers
           times the rows of U; its entries between the columns rank and
           col are zero already (they were zeroed when these columns were
           processed) */
        if (lo[rank] < hi[rank] && col + 1 < n)
        {
            slong l = (lo[rank] < PLE_CROUT_MIN_SKIP) ? 0 : lo[rank];
            slong h = hi[rank];
            nmod_mat_t Uw;

            nmod_mat_window_init(Uw, A, l, col + 1, h, n);
            nmod_mat_nmod_vec_mul(tmp, prow + l, h - l, Uw);
            nmod_mat_window_clear(Uw);
            _nmod_vec_sub(prow + col + 1, prow + col + 1, tmp, n - col - 1, mod);
        }

        inv = nmod_inv(prow[col], mod);

        /* the multipliers, in column rank */
        for (i = rank + 1; i < m; i++)
        {
            ulong * row = ROW(i);
            ulong e = row[col];

            row[col] = 0;
            row[rank] = (inv == 1) ? e : nmod_mul(e, inv, mod);
            if (e != 0)
            {
                if (lo[i] >= hi[i])
                    lo[i] = rank;
                hi[i] = rank + 1;
            }
        }

        Q[rank] = col;
        rank++;
    }

#undef ROW

    _nmod_mat_ple_rows_finish(P, Q, A, orig, rank);

    TMP_END;

    return rank;
}

slong
nmod_mat_ple_crout(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
{
    return _nmod_mat_ple_wrap(P, pivots_nonpivots, A, flags, _nmod_mat_ple_crout);
}
