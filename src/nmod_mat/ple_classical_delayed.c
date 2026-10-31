/*
    Copyright (C) 2021 Fredrik Johansson
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
#include "nmod_mat/impl.h"
#include "nmod_mat/delayed_impl.h"

/*
    The elimination of ple_classical.c with the reductions of
    lu_classical_delayed.c: the updates of an entry are accumulated
    unreduced, in A itself (one word) or in a buffer of two or three words
    per entry, and an entry is reduced when it is about to be used: when
    its column is searched for a pivot (rows that are not pivot rows yet)
    and when its row becomes the pivot row. Before the first pivot is
    found, no update has been made and nothing needs to be reduced. Rows
    with a zero multiplier are not updated, so that the multipliers added
    are in [1, p - 1] and every unreduced entry is at most s (p - 1)^2,
    p the modulus and s = min(m, n).

    As in nmod_mat_lu_classical_delayed, a new pivot row is swapped with
    the row of A (and of the buffer) at the position of the rank, so that
    the eliminations run through consecutive rows. The rows that are not
    pivot rows are thus stored in any order; the pivot is the one with a
    nonzero entry of smallest index in the input (orig), found while the
    column is reduced, and their order is restored at the end (see
    _nmod_mat_ple_rows_finish in impl.h).
*/

/* the constants of the reductions, for words = 1, 2, 3 words per entry */
typedef struct
{
    nmod_t mod;
    ulong npre, npre2;      /* one word: Barrett, Lemire */
    int lemire;             /* one word: unreduced entries below 2^(FLINT_BITS/2) */
    int lemire_mul;         /* one word: products of reduced entries likewise */
    int halflimb;           /* two words: n <= 2^(FLINT_BITS/2) */
    int high_reduced;       /* two words: high words below n */
    ulong qhi, qlo;         /* two words: floor(2^(2 FLINT_BITS) / n) */
    int fullword_limited;   /* three words: n close to 2^(FLINT_BITS-1) */
    ulong alpha1, alpha2;   /* three words: 2^FLINT_BITS, 2^(2 FLINT_BITS) mod n */
}
_ple_red_ctx;

/* the entry x (words words) reduced */
FLINT_FORCE_INLINE ulong
_ple_red(const ulong * x, int words, const _ple_red_ctx * C)
{
    ulong r;

    if (words == 1)
    {
        if (C->lemire)
            r = n_mod_lemire(x[0], C->mod.n, C->npre2);
        else
            r = n_mod_barrett(x[0], C->mod.n, C->npre);
    }
    else if (words == 2)
    {
        if (C->high_reduced)
            NMOD_RED2(r, x[1], x[0], C->mod);
        else
            r = n_ll_rem_l_nonfullword(x[1], x[0], C->mod.n, C->qhi, C->qlo);
    }
    else
    {
        if (C->fullword_limited)
            r = n_lll_rem_l_fullword_limited(x[2], x[1], x[0], C->mod,
                                             C->alpha2, C->alpha1);
        else
            r = n_lll_rem_l(x[2], x[1], x[0], C->mod, C->alpha2, C->alpha1);
    }

    return r;
}

/* x[0..len) += f * v[0..len), unreduced */
FLINT_FORCE_INLINE void
_ple_addmul(ulong * x, const ulong * v, slong len, ulong f, int words,
            const _ple_red_ctx * C)
{
    if (words == 1)
        _nmod_vec_nored_scalar_addmul_halflimb(x, v, len, f);
    else if (words == 2)
    {
        if (C->halflimb)
            _nmod_vec_nored_ll_scalar_addmul_halflimb(x, v, len, f);
        else
            _nmod_vec_nored_ll_scalar_addmul(x, v, len, f);
    }
    else
        _nmod_vec_nored_lll_scalar_addmul(x, v, len, f);
}

/*
    The elimination with words words per unreduced entry: in A itself for
    one word (b is NULL), otherwise in b, the entry (i, j) of A having its
    words at b + words (i n + j).
*/
FLINT_FORCE_INLINE slong
_ple_delayed(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots,
             int words, ulong * b, const _ple_red_ctx * C)
{
    const slong m = A->r, n = A->c, stride = A->stride;
    ulong * const aa = A->entries;
    const nmod_t mod = C->mod;
    slong rank = 0, nonpivots = 0, col, i, j, pr;
    slong * orig;
    TMP_INIT;

    /* local copies of the entries and the stride: the stores to the
       entries cannot change them */
#define ROW(i) (aa + (i) * stride)
#define UNRED(i, j) ((words == 1) ? ROW(i) + (j) : b + words * ((i) * n + (j)))

    TMP_START;
    orig = TMP_ALLOC((m + 1) * sizeof(slong));

    for (i = 0; i < m; i++)
        orig[i] = i;

    for (col = 0; col < n && rank < m; col++)
    {
        ulong * prow;
        ulong d;

        /* the pivot: among the rows that are not pivot rows, the one with
           a nonzero entry in this column of smallest index in the input;
           before the first pivot, the rows are in the order of the input
           and need no reduction */
        pr = -1;
        if (rank > 0)
        {
            slong best = m;

            for (i = rank; i < m; i++)
            {
                ulong x = _ple_red(UNRED(i, col), words, C);
                slong o = orig[i];
                int better = (x != 0) & (o < best);

                ROW(i)[col] = x;
                /* without branches: zero entries can be frequent and
                   unpredictable */
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
        }
        P[rank] = orig[rank];

        /* the pivot row, reduced; with a buffer, the other row moves to
           row pr of the buffer (its entries up to column col are no longer
           used) */
        prow = ROW(rank);
        if (words == 1)
        {
            if (rank > 0)
                for (j = col + 1; j < n; j++)
                    prow[j] = _ple_red(prow + j, 1, C);
        }
        else if (pr != rank)
        {
            for (j = col + 1; j < n; j++)
            {
                ulong * x = UNRED(pr, j);
                ulong * y = UNRED(rank, j);

                if (rank > 0)
                    prow[j] = _ple_red(x, words, C);
                x[0] = y[0];
                x[1] = y[1];
                if (words == 3)
                    x[2] = y[2];
            }
        }
        else if (rank > 0)
        {
            for (j = col + 1; j < n; j++)
                prow[j] = _ple_red(UNRED(rank, j), words, C);
        }

        d = nmod_inv(prow[col], mod);

        for (i = rank + 1; i < m; i++)
        {
            ulong * row = ROW(i);
            ulong e = row[col];

            if (e == 0)
                continue;

            if (words == 1)
            {
                if (d != 1)
                    e = C->lemire_mul ? n_mod_lemire(d * e, mod.n, C->npre2)
                                      : n_mod_barrett(d * e, mod.n, C->npre);
            }
            else
                e = nmod_mul(d, e, mod);

            _ple_addmul(UNRED(i, col + 1), prow + col + 1, n - col - 1,
                        mod.n - e, words, C);
            row[col] = 0;
            row[rank] = e;
        }

        Q[rank] = col;
        rank++;
    }

#undef ROW
#undef UNRED

    _nmod_mat_ple_rows_finish(P, Q, A, orig, rank);

    TMP_END;

    return rank;
}

static slong
_ple_delayed_1(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    _ple_red_ctx C;
    ulong bound;

    memset(&C, 0, sizeof(C));
    C.mod = A->mod;
    bound = (C.mod.n - 1) * (C.mod.n - 1) * FLINT_MIN(A->r, A->c);
    C.lemire = bound < (UWORD(1) << (FLINT_BITS / 2));
    C.lemire_mul = NMOD_BITS(C.mod) <= FLINT_BITS / 4;
    C.npre = n_barrett_precomp(C.mod.n);
    C.npre2 = n_lemire_precomp(C.mod.n);

    return _ple_delayed(P, Q, A, max_nonpivots, 1, NULL, &C);
}

static slong
_ple_delayed_23(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots,
                int words)
{
    const slong m = A->r, n = A->c;
    _ple_red_ctx C;
    ulong * b;
    slong i, j, rank;
    TMP_INIT;

    memset(&C, 0, sizeof(C));
    C.mod = A->mod;

    if (words == 2)
    {
        /* only reachable with min(m, n) <= 3 for p >= 2^(FLINT_BITS - 1) */
        if (C.mod.norm == 0)
            return _nmod_mat_ple_classical(P, Q, A, max_nonpivots);

        C.halflimb = NMOD_IS_HALFLIMB(C.mod);
        C.high_reduced = 0;
        if (!C.halflimb)
        {
            ulong hi, lo;
            umul_ppmm(hi, lo, C.mod.n - 1, C.mod.n - 1);
            mullo_2x1(&hi, &lo, hi, lo, FLINT_MIN(m, n));
            C.high_reduced = (hi < C.mod.n);
        }
        C.qhi = C.qlo = 0;
        if (!C.high_reduced)
            n_ll_rem_l_precomp(&C.qhi, &C.qlo, C.mod.n);
    }
    else
    {
        /* moduli close to 2^(FLINT_BITS - 1) */
        C.fullword_limited = (C.mod.norm == 0) && C.mod.n <
                (UWORD(1) << (FLINT_BITS - 1)) + (UWORD(1) << (FLINT_BITS / 2 - 2));

        if (C.fullword_limited)
        {
            C.alpha1 = -C.mod.n;                 /* 2^FLINT_BITS */
            C.alpha2 = 4 * C.alpha1 * C.alpha1;  /* 2^(2 FLINT_BITS) */
        }
        else
        {
            C.alpha1 = nmod_set_ui(UWORD(1) << (FLINT_BITS - 1), C.mod);
            C.alpha1 = nmod_add(C.alpha1, C.alpha1, C.mod);  /* 2^FLINT_BITS */
            C.alpha2 = nmod_mul(C.alpha1, C.alpha1, C.mod);  /* 2^(2 FLINT_BITS) */
        }
    }

    TMP_START;
    b = TMP_ALLOC(words * m * n * sizeof(ulong));

    for (i = 0; i < m; i++)
    {
        const ulong * row = nmod_mat_entry_ptr(A, i, 0);
        for (j = 0; j < n; j++)
        {
            b[words * (i * n + j)] = row[j];
            b[words * (i * n + j) + 1] = 0;
            if (words == 3)
                b[words * (i * n + j) + 2] = 0;
        }
    }

    if (words == 2)
        rank = _ple_delayed(P, Q, A, max_nonpivots, 2, b, &C);
    else
        rank = _ple_delayed(P, Q, A, max_nonpivots, 3, b, &C);

    TMP_END;

    return rank;
}

slong
_nmod_mat_ple_classical_delayed(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    int nlimbs;

    /* for the modulus 1, all entries are zero and the reductions above
       cannot be precomputed */
    if (A->mod.n == 1)
        return _nmod_mat_ple_classical(P, Q, A, max_nonpivots);

    /* number of limbs of the bound s (p - 1)^2 on the unreduced entries */
    nlimbs = _nmod_vec_dot_bound_limbs(FLINT_MIN(A->r, A->c), A->mod);

    if (nlimbs <= 1)
        return _ple_delayed_1(P, Q, A, max_nonpivots);
    else if (nlimbs == 2)
        return _ple_delayed_23(P, Q, A, max_nonpivots, 2);
    else
        return _ple_delayed_23(P, Q, A, max_nonpivots, 3);
}

slong
nmod_mat_ple_classical_delayed(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
{
    return _nmod_mat_ple_wrap(P, pivots_nonpivots, A, flags,
                              _nmod_mat_ple_classical_delayed);
}
