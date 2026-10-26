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
#include "nmod_mat/impl.h"

/*
    Rows of A (all its columns) permuted in place so that the new row i is
    the old row perm[i], one cycle of perm at a time with one temporary
    row; rows fixed by perm are not touched.
*/
void
_nmod_mat_ple_permute_rows(nmod_mat_t A, const slong * perm)
{
    slong m = A->r, n = A->c, i, j, start;
    ulong * tmp;
    char * done;
    TMP_INIT;

    if (n <= 0)
        return;

    for (start = 0; start < m && perm[start] == start; start++)
        ;

    if (start >= m)
        return;

    TMP_START;
    tmp = TMP_ALLOC(n * sizeof(ulong));
    done = TMP_ALLOC(m * sizeof(char));
    memset(done, 0, m);

    for (i = start; i < m; i++)
    {
        if (done[i] || perm[i] == i)
            continue;

        memcpy(tmp, nmod_mat_entry_ptr(A, i, 0), n * sizeof(ulong));

        for (j = i; perm[j] != i; j = perm[j])
        {
            memcpy(nmod_mat_entry_ptr(A, j, 0),
                   nmod_mat_entry_ptr(A, perm[j], 0), n * sizeof(ulong));
            done[j] = 1;
        }

        memcpy(nmod_mat_entry_ptr(A, j, 0), tmp, n * sizeof(ulong));
        done[j] = 1;
    }

    TMP_END;
}

/*
    End of the base cases, which work with the rows in place and the
    logical order in P: moves row P[i] to row i, and writes the columns
    without pivot in increasing order to Q[rank..n), after the pivot
    columns Q[0..rank) (increasing).
*/
void
_nmod_mat_ple_finish(slong * P, slong * Q, nmod_mat_t A, slong rank)
{
    slong n = A->c, j, k, l;

    for (j = k = 0, l = rank; j < n; j++)
    {
        if (k < rank && Q[k] == j)
            k++;
        else
            Q[l++] = j;
    }

    _nmod_mat_ple_permute_rows(A, P);
}

void
_nmod_mat_ple_rows_finish(slong * P, slong * Q, nmod_mat_t A, slong * orig,
                          slong rank)
{
    slong m = A->r, n = A->c, i, j, k;

    if (rank < m)
    {
        slong * pos;
        nmod_mat_t W;
        TMP_INIT;

        TMP_START;
        pos = TMP_ALLOC(m * sizeof(slong));

        /* P[rank..m): the rows that are not pivot rows, increasing */
        for (i = 0; i < m; i++)
            pos[i] = -1;
        for (k = 0; k < rank; k++)
            pos[P[k]] = k;
        for (i = 0, k = rank; i < m; i++)
            if (pos[i] < 0)
                P[k++] = i;

        /* their rows in A, in this order */
        for (i = rank; i < m; i++)
            pos[orig[i]] = i;
        for (k = rank; k < m; k++)
            orig[k - rank] = pos[P[k]] - rank;
        nmod_mat_window_init(W, A, rank, 0, m, n);
        _nmod_mat_ple_permute_rows(W, orig);
        nmod_mat_window_clear(W);

        TMP_END;
    }

    for (j = k = 0, i = rank; j < n; j++)
    {
        if (k < rank && Q[k] == j)
            k++;
        else
            Q[i++] = j;
    }
}

/*
    The choice of the algorithm, tuned on uniformly random matrices (one
    thread, no external BLAS) on Zen 4 and Ice Lake SP (AVX512-IFMA),
    Cascade Lake (AVX-512 without IFMA) and Apple M4 (NEON); the profile
    program p-ple gives the timings to tune it on other machines.

    The crossovers depend on the SIMD kernels compiled in, hence two classes:
    with AVX512-IFMA (NMOD_MAT_HAVE_MUL_U52), whose u52 dot products and
    vector-matrix products make the Crout algorithm faster up to 2^52, and the
    others (AVX2 or AVX-512 without IFMA, NEON, generic), with two values for
    30 to 32 bits that differ between x86 (AVX2, NMOD_MAT_HAVE_VECMAT_U64) and
    the other targets. 

    TODO AVX2 measurements (the second class was measured with AVX-512 and
    NEON only: for AVX2 alone the values are extrapolated)

    - min(m, n) <= 3, or <= 7 when the delayed algorithm would need two or
      three words per entry: the classical algorithm, as in nmod_mat_lu,
      except for min(m, n) = 2 or 3 with one word per entry and a matrix
      more than 4 times wider than tall, or more than 2 times (IFMA) or 8
      times (others) taller than wide (not for p = 2). Without IFMA, also
      for n <= 8 and n / 2 <= m <= 2 n above 20 bits (at least n rows up to
      24 bits, at most 8 rows up to 32 bits), where the classical algorithm
      is inversion-free (see _nmod_mat_ple_classical).
    - The modulus 2 on x86 with AVX2, for m x n with n >= 4 m below the
      recursion: the classical algorithm too (its row operations are
      additions; 1.2-1.75x faster than the delayed algorithm on the three
      x86 machines, 1.15-2.2x slower on the M4).
    - The recursive algorithm when n >= N, m >= M and m n >= 2^e, with
      (N, M, e) from ple_rec_wanted, by class (IFMA, NEON, others) and
      modulus size. The recursion splits the columns: for a given m n it pays
      off earlier for tall matrices than for wide ones, hence the bounds on
      n and m. Fitted on the crossovers between the base case and one level
      of recursion, on Zen 4 and Ice Lake SP (IFMA), Apple M4 (NEON) and
      Cascade Lake (others), then checked with the whole recursion, also on
      matrices of rank n / 2, on which the recursion pays off later (hence
      the higher bounds up to 29 bits). Without IFMA, recursion is used even
      above 2^52 for large matrices (m n >= 2^19).
    - Below, the Crout algorithm when the delayed algorithm would need two
      or three words per entry, nmod_mat_nmod_vec_mul has a SIMD path for
      the modulus (not above 2^50 on NEON, where Crout was up to 1.5x
      slower), and min(m, n) >= S or (n >= 2 m and n >= N0), with (S, N0)
      from ple_crout_min (IFMA / other): (32, 16) / (64, 256) up to 32 bits
      ((128, 256) without AVX2),
      (32, 16) / (128, 64) from 33 to 50 bits, (32, 16) / (128, 128) for 51
      and 52 bits, (64, 128) / (128, 128) from 53 to 62 bits, (64, 128) /
      (64, 128) above. Otherwise the delayed algorithm.

    The classical and delayed algorithms skip the rows whose multiplier is
    zero, the Crout algorithm only part of them (at the ends of each row)
    and the recursion none: on matrices with many zero multipliers
    (structured or sparse matrices, non-generic rank profiles), the delayed
    algorithm can be faster than the others beyond these crossovers. On the
    L R U matrices of p-ple (ranks), whose L has many zero multipliers
    scattered in each row, the choice above is 1.09-1.17x slower than the
    best on average (1.08-1.15x before this tuning, which did not use
    Crout as much), and up to 2-2.6x for some low ranks.
*/

/* whether the classical algorithm is used for small m x n */
static int
ple_classical_wanted(slong m, slong n, flint_bitcnt_t bits, int nlimbs,
                     ulong p)
{
    slong s = FLINT_MIN(m, n);
#if NMOD_MAT_HAVE_MUL_U52
    const slong tall = 2;
#else
    const slong tall = 8;
#endif

    if (s <= 3)
        return s <= 1 || p == 2 || nlimbs >= 2
            || (m <= tall * n && n <= 4 * m);

    if (s <= 7 && nlimbs >= 2)
        return 1;

#if !NMOD_MAT_HAVE_MUL_U52
    if (n <= 8 && m <= 2 * n && 2 * m >= n && (m >= n || bits > 24)
            && (m <= 8 || bits > 32) && bits > 20)
        return 1;
#else
    (void) bits;
#endif

    return 0;
}

/* whether the recursive algorithm is used for m x n: n >= N, m >= M and
   m n >= 2^e, with (N, M, e) by class and modulus size */
static int
ple_rec_wanted(slong m, slong n, flint_bitcnt_t bits)
{
    slong N, M;
    int e;

#if NMOD_MAT_HAVE_MUL_U52
    if (bits <= 2)       N = 128, M = 128, e = 14;
    else if (bits <= 29) N = 128, M = 128, e = 15;
    else if (bits <= 30) N = 64, M = 64, e = 0;
    else if (bits <= 32) N = 64, M = 128, e = 0;
    else if (bits <= 52) N = 16, M = 128, e = 12;
    else                 return 0;
#elif (defined(__ARM_NEON) && defined(__aarch64__)) || defined(_M_ARM64)
    if (bits <= 29)      N = 128, M = 128, e = 15;
    else if (bits <= 30) N = 64, M = 64, e = 0;
    else if (bits <= 32) N = 96, M = 96, e = 14;
    else if (bits <= 52) N = 96, M = 96, e = 16;
    else                 N = 96, M = 96, e = 19;
#else
    if (bits <= 2)       N = 128, M = 128, e = 14;
    else if (bits <= 29) N = 128, M = 128, e = 15;
    else if (bits <= 30) N = 32, M = 32, e = 13;
    else if (bits <= 32) N = 16, M = 128, e = 11;
    else if (bits <= 52) N = 16, M = 96, e = 12;
    else                 N = 32, M = 32, e = 19;
#endif

    return n >= N && m >= M && (ulong) m * (ulong) n >= (UWORD(1) << e);
}

/* (S, N0) of the condition for the Crout algorithm */
static void
ple_crout_min(slong * S, slong * N0, flint_bitcnt_t bits)
{
#if NMOD_MAT_HAVE_MUL_U52
    if (bits <= 52)
        *S = 32, *N0 = 16;
    else
        *S = 64, *N0 = 128;
#else
    if (bits <= 32)
        *S = NMOD_MAT_HAVE_VECMAT_U64 ? 64 : 128, *N0 = 256;
    else if (bits <= 50)
        *S = 128, *N0 = 64;
    else if (bits <= 62)
        *S = 128, *N0 = 128;
    else
        *S = 64, *N0 = 128;
#endif
}

slong
_nmod_mat_ple(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    slong m = A->r, n = A->c, s = FLINT_MIN(m, n), S, N0;
    flint_bitcnt_t bits = NMOD_BITS(A->mod);
    int nlimbs = _nmod_vec_dot_bound_limbs(s, A->mod);

    if (ple_classical_wanted(m, n, bits, nlimbs, A->mod.n))
        return _nmod_mat_ple_classical(P, Q, A, max_nonpivots);

    if (ple_rec_wanted(m, n, bits))
        return _nmod_mat_ple_recursive(P, Q, A, max_nonpivots);

#if NMOD_MAT_HAVE_VECMAT_U64
    if (A->mod.n == 2 && n >= 4 * m)
        return _nmod_mat_ple_classical(P, Q, A, max_nonpivots);
#endif

    ple_crout_min(&S, &N0, bits);
    if (nlimbs >= 2 && NMOD_MAT_NMOD_VEC_MUL_IS_SIMD(A->mod.n)
            && (s >= S || (n >= 2 * m && n >= N0)))
        return _nmod_mat_ple_crout(P, Q, A, max_nonpivots);

    return _nmod_mat_ple_classical_delayed(P, Q, A, max_nonpivots);
}

slong
_nmod_mat_ple_wrap(slong * P, slong * pivots_nonpivots, nmod_mat_t A,
        int flags, slong (* ple)(slong *, slong *, nmod_mat_t, slong))
{
    slong m = A->r, n = A->c, rank, max_nonpivots;
    slong * Q = pivots_nonpivots;

    /* with the rank check, stop as soon as more than n - min(m, n)
       columns without pivot are met: the rank is then less than min(m, n) */
    max_nonpivots = (flags & NMOD_MAT_PLE_RANK_CHECK) ? n - FLINT_MIN(m, n) : -1;

    if (Q == NULL)
        Q = flint_malloc((n + 1) * sizeof(slong));

    rank = ple(P, Q, A, max_nonpivots);

    if (pivots_nonpivots == NULL)
        flint_free(Q);

    return FLINT_MAX(rank, 0);
}

slong
nmod_mat_ple(slong * P, slong * pivots_nonpivots, nmod_mat_t A, int flags)
{
    return _nmod_mat_ple_wrap(P, pivots_nonpivots, A, flags, _nmod_mat_ple);
}
