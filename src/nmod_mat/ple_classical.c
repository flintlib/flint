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
    Whether the inversion-free variant below is used, which pays when there
    are few entries to update per inversion saved: moduli of at least 16
    bits (below, the inversions are cheap), rows of at most PLE_NOINV_MAX
    entries (for longer rows, the vectorized updates of the classical
    algorithm win), at least n rows up to 24 bits and at least n / 2 rows
    up to 32 bits (with fewer rows, fewer row updates per inversion saved);
    for full-word moduli, whose updates cost two reductions, at most 16 rows.
*/
#define PLE_NOINV_MAX 8

static int
ple_noinv_wanted(slong m, slong n, nmod_t mod)
{
    const flint_bitcnt_t bits = NMOD_BITS(mod);

    return n <= PLE_NOINV_MAX && bits >= 16
        && (m >= n || bits > 24) && (2 * m >= n || bits > 32)
        && (mod.norm > 0 || m <= 16);
}

/*
    The same elimination without inversions until the end, for small
    matrices, where inverting each pivot costs about as much as the updates
    of the entries. Each row that is not a pivot row yet holds sigma times
    its value in the algorithm below, with the same factor sigma for all of
    them: eliminating with the pivot row (entry d in the pivot column) does
    R_i <- d R_i - e R_pivot (or R_i <- d R_i if e = 0), and sigma <- d sigma.
    The multipliers e and the pivot rows are thus scaled, by the factor
    sigma of their step; one inversion of the last factor gives the inverses
    of the others and of the pivots (d = sigma' / sigma), which normalize
    them at the end. For a composite modulus, a pivot that is not
    invertible makes the last factor not invertible, and the exception is
    raised at the end.
*/

/* row[0..len) <- d row - e prow, with ne = n - e */
FLINT_FORCE_INLINE void
_ple_noinv_update(ulong * row, const ulong * prow, slong len, ulong d,
                  ulong ne, nmod_t mod)
{
    slong j;

    if (mod.norm > 0)
    {
        /* the sum of the two products is below n 2^FLINT_BITS */
        for (j = 0; j < len; j++)
        {
            ulong h1, l1, h2, l2;
            umul_ppmm(h1, l1, d, row[j]);
            umul_ppmm(h2, l2, ne, prow[j]);
            add_ssaaaa(h1, l1, h1, l1, h2, l2);
            NMOD_RED2(row[j], h1, l1, mod);
        }
    }
    else
    {
        for (j = 0; j < len; j++)
            row[j] = nmod_add(nmod_mul(d, row[j], mod),
                              nmod_mul(ne, prow[j], mod), mod);
    }
}

/* out of line, so that it does not weigh on the other function */
FLINT_STATIC_NOINLINE slong
_ple_classical_noinv(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    const slong m = A->r, n = A->c;
    const nmod_t mod = A->mod;
    slong rank = 0, nonpivots = 0, col, i, j, r;
    ulong sig[PLE_NOINV_MAX + 1], D[PLE_NOINV_MAX], isig[PLE_NOINV_MAX + 1];

    FLINT_ASSERT(n <= PLE_NOINV_MAX);

    for (i = 0; i < m; i++)
        P[i] = i;

    sig[0] = 1;

    for (col = 0; col < n && rank < m; col++)
    {
        ulong * prow;
        ulong d;
        slong len;

        for (r = rank; r < m && nmod_mat_entry(A, P[r], col) == 0; r++)
            ;

        if (r == m)
        {
            nonpivots++;
            if (max_nonpivots >= 0 && nonpivots > max_nonpivots)
            {
                /* as the other algorithms, raise the exception if a pivot
                   met so far is not invertible */
                if (rank > 0)
                    nmod_inv(sig[rank], mod);
                return -1;
            }
            continue;
        }

        _nmod_mat_ple_rotate(P, rank, r);

        prow = nmod_mat_entry_ptr(A, P[rank], 0);
        d = prow[col];
        len = n - col - 1;

        for (i = rank + 1; i < m; i++)
        {
            ulong * row = nmod_mat_entry_ptr(A, P[i], 0);
            ulong e = row[col];

            if (e == 0)
            {
                /* row[col] and row[rank] are already zero */
                if (d != 1)
                    for (j = col + 1; j < n; j++)
                        row[j] = nmod_mul(d, row[j], mod);
                continue;
            }

            _ple_noinv_update(row + col + 1, prow + col + 1, len, d,
                              mod.n - e, mod);
            row[col] = 0;
            row[rank] = e;
        }

        D[rank] = d;
        sig[rank + 1] = nmod_mul(sig[rank], d, mod);
        Q[rank] = col;
        rank++;
    }

    if (rank > 0)
    {
        /* isig[t] = 1 / sigma_t, the multipliers of column t are divided
           by D[t] = sigma_{t+1} / sigma_t, the pivot row t by sigma_t */
        isig[rank] = nmod_inv(sig[rank], mod);
        for (r = rank - 1; r >= 1; r--)
            isig[r] = nmod_mul(isig[r + 1], D[r], mod);
        isig[0] = 1;
        for (r = 0; r < rank; r++)
            D[r] = nmod_mul(sig[r], isig[r + 1], mod);

        for (i = 0; i < m; i++)
        {
            ulong * row = nmod_mat_entry_ptr(A, P[i], 0);
            slong t = FLINT_MIN(i, rank);

            for (j = 0; j < t; j++)
                row[j] = nmod_mul(row[j], D[j], mod);

            if (i > 0 && i < rank)
                for (j = Q[i]; j < n; j++)
                    row[j] = nmod_mul(row[j], isig[i], mod);
        }
    }

    _nmod_mat_ple_finish(P, Q, A, rank);

    return rank;
}

/*
    Gaussian elimination column by column, with a reduction after each
    update. The rows stay in place: the row in logical position i is row
    P[i] of A, and bringing the pivot row to position rank is an upward rotation
    of P[rank..r] (the rows passed over keep their relative order). Rows with
    a zero multiplier are not updated.
*/
slong
_nmod_mat_ple_classical(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    const slong m = A->r, n = A->c;
    const nmod_t mod = A->mod;
    slong rank = 0, nonpivots = 0, col, i, r;

    if (ple_noinv_wanted(m, n, mod))
        return _ple_classical_noinv(P, Q, A, max_nonpivots);

    for (i = 0; i < m; i++)
        P[i] = i;

    for (col = 0; col < n && rank < m; col++)
    {
        ulong * prow;
        ulong d;
        slong len;

        /* pivot: first row in logical order with a nonzero entry */
        for (r = rank; r < m && nmod_mat_entry(A, P[r], col) == 0; r++)
            ;

        if (r == m)
        {
            nonpivots++;
            if (max_nonpivots >= 0 && nonpivots > max_nonpivots)
                return -1;
            continue;
        }

        _nmod_mat_ple_rotate(P, rank, r);

        prow = nmod_mat_entry_ptr(A, P[rank], 0);
        d = nmod_inv(prow[col], mod);
        len = n - col - 1;

        for (i = rank + 1; i < m; i++)
        {
            ulong * row = nmod_mat_entry_ptr(A, P[i], 0);
            ulong e = row[col];

            /* nothing to do: row[col] and row[rank] are already zero */
            if (e == 0)
                continue;

            e = nmod_mul(e, d, mod);
            if (len > 0)
                _nmod_vec_scalar_addmul_nmod(row + col + 1, prow + col + 1,
                                             len, nmod_neg(e, mod), mod);
            row[col] = 0;
            row[rank] = e;
        }

        Q[rank] = col;
        rank++;
    }

    _nmod_mat_ple_finish(P, Q, A, rank);

    return rank;
}

slong
nmod_mat_ple_classical(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
{
    return _nmod_mat_ple_wrap(P, pivots_nonpivots, A, flags,
                              _nmod_mat_ple_classical);
}
