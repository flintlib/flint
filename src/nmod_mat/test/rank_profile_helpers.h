/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* shared helpers of the tests involving rank profiles */

#ifndef NMOD_MAT_TEST_RANK_PROFILE_HELPERS_H
#define NMOD_MAT_TEST_RANK_PROFILE_HELPERS_H

#include "ulong_extras.h"
#include "perm.h"
#include "nmod.h"
#include "nmod_mat.h"

/* comparison of two slong, for flint_sort */
static int
_rp_slong_cmp(const void * a, const void * b, void * FLINT_UNUSED(data))
{
    slong x = *(const slong *) a, y = *(const slong *) b;

    return (x > y) - (x < y);
}

/*
    Random positions (mrp_r[k], mrp_c[k]), k < rank, in distinct rows and
    columns of an m x n matrix, sorted by increasing column. Several
    shapes: random; on the diagonal (generic rank profile); rows
    increasing with the columns (the pivot rows come in increasing
    order); rows decreasing (each new pivot lies above all the previous
    ones, so that each pivot row passes over all the rows above it that
    are not pivot rows yet).
*/
static void
_rp_rand_positions(slong * mrp_r, slong * mrp_c, slong m, slong n,
                                    slong rank, flint_rand_t state)
{
    slong k;
    slong * pr = _perm_init(m);
    slong * pc = _perm_init(n);
    int shape = n_randint(state, 4);

    _perm_randtest(pr, m, state);
    _perm_randtest(pc, n, state);

    /* rank random columns in increasing order, each with a random row */
    flint_sort(pc, rank, sizeof(slong), _rp_slong_cmp, NULL);
    for (k = 0; k < rank; k++)
    {
        mrp_r[k] = pr[k];
        mrp_c[k] = pc[k];
    }

    if (shape == 1)
    {
        for (k = 0; k < rank; k++)
            mrp_r[k] = mrp_c[k] = k;
    }
    else if (shape >= 2)
    {
        flint_sort(mrp_r, rank, sizeof(slong), _rp_slong_cmp, NULL);
        if (shape == 3)
            for (k = 0; k < rank / 2; k++)
                FLINT_SWAP(slong, mrp_r[k], mrp_r[rank - 1 - k]);
    }

    _perm_clear(pr);
    _perm_clear(pc);
}

/*
    Rank profile matrix of A by its definition, for small matrices and a
    prime modulus: there is a pivot at (i, j) if and only if
    rho(i + 1, j + 1) - rho(i, j + 1) - rho(i + 1, j) + rho(i, j) = 1,
    where rho(a, b) is the rank of the leading a x b submatrix. Writes the
    pivots sorted by increasing column to mrp_r and mrp_c, returns the
    rank.
*/
static slong
_rp_rank_profile_naive(slong * mrp_r, slong * mrp_c, const nmod_mat_t A)
{
    slong m = A->r, n = A->c, i, j, r;
    slong * rho = flint_calloc((m + 1) * (n + 1), sizeof(slong));
    nmod_mat_t W;

    for (i = 1; i <= m; i++)
        for (j = 1; j <= n; j++)
        {
            nmod_mat_window_init(W, A, 0, 0, i, j);
            rho[i * (n + 1) + j] = nmod_mat_rank(W);
            nmod_mat_window_clear(W);
        }

    r = 0;
    for (j = 0; j < n; j++)
        for (i = 0; i < m; i++)
            if (rho[(i + 1) * (n + 1) + j + 1] - rho[i * (n + 1) + j + 1]
                    - rho[(i + 1) * (n + 1) + j] + rho[i * (n + 1) + j] == 1)
            {
                mrp_r[r] = i;
                mrp_c[r] = j;
                r++;
            }

    flint_free(rho);
    return r;
}

/*
    Checks that (P, Q, LU, rank) is the canonical PLE decomposition of A
    (see nmod_mat_ple): LU in the storage of nmod_mat_lu, row i of L U
    equal to row P[i] of A, Q = pivot columns of U (increasing) followed
    by the other columns (increasing), P[rank..m) increasing, and rank
    profile matrix {(P[k], Q[k]), k < rank}. If mrp_r is not NULL, the
    rank profile matrix of A is given by mrp_r, mrp_c (sorted by column) and its
    rank by exp_rank; otherwise, if check_rpm is nonzero, it is computed by
    _rp_rank_profile_naive (small matrices). Returns 0 if correct, a
    positive code otherwise.
*/
static int
_rp_check_ple(const slong * P, const slong * Q, const nmod_mat_t LU,
              slong rank, const nmod_mat_t A, const slong * mrp_r,
              const slong * mrp_c, slong exp_rank, int check_rpm)
{
    slong m = A->r, n = A->c, i, j, k;
    int res = 0;
    char * seen;
    nmod_mat_t L, U, B;

    if (rank < 0 || rank > FLINT_MIN(m, n))
        return 1;

    /* P and Q are permutations, with the monotonicity conditions */
    seen = flint_calloc(m + n + 1, sizeof(char));
    for (i = 0; i < m && res == 0; i++)
    {
        if (P[i] < 0 || P[i] >= m || seen[P[i]])
            res = 2;
        else
            seen[P[i]] = 1;
        if (i > rank && P[i] <= P[i - 1])
            res = 3;
    }
    for (j = 0; j < n && res == 0; j++)
    {
        if (Q[j] < 0 || Q[j] >= n || seen[m + Q[j]])
            res = 4;
        else
            seen[m + Q[j]] = 1;
        if (j > 0 && j != rank && Q[j] <= Q[j - 1])
            res = 5;
    }
    flint_free(seen);
    if (res != 0)
        return res;

    /* shape: U in echelon form with pivots in the columns Q[0..rank),
       rows rank..m-1 zero beyond column rank */
    for (k = 0; k < rank; k++)
    {
        for (j = k; j < Q[k]; j++)
            if (nmod_mat_entry(LU, k, j) != 0)
                return 6;
        if (nmod_mat_entry(LU, k, Q[k]) == 0)
            return 7;
    }
    for (i = rank; i < m; i++)
        for (j = rank; j < n; j++)
            if (nmod_mat_entry(LU, i, j) != 0)
                return 8;

    /* L U = rows of A permuted by P */
    nmod_mat_init(L, m, rank, A->mod.n);
    nmod_mat_init(U, rank, n, A->mod.n);
    nmod_mat_init(B, m, n, A->mod.n);
    for (i = 0; i < m; i++)
        for (k = 0; k < FLINT_MIN(i + 1, rank); k++)
            nmod_mat_entry(L, i, k) = (i == k) ? 1 : nmod_mat_entry(LU, i, k);
    for (k = 0; k < rank; k++)
        for (j = k; j < n; j++)
            nmod_mat_entry(U, k, j) = nmod_mat_entry(LU, k, j);
    nmod_mat_mul(B, L, U);
    for (i = 0; i < m && res == 0; i++)
        for (j = 0; j < n; j++)
            if (nmod_mat_entry(B, i, j) != nmod_mat_entry(A, P[i], j))
            {
                res = 9;
                break;
            }
    nmod_mat_clear(L);
    nmod_mat_clear(U);
    nmod_mat_clear(B);
    if (res != 0)
        return res;

    /* rank profile matrix */
    if (mrp_r != NULL || check_rpm)
    {
        slong r2 = exp_rank;
        slong * mrp_r2 = (slong *) mrp_r;
        slong * mrp_c2 = (slong *) mrp_c;

        if (mrp_r == NULL)
        {
            mrp_r2 = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
            mrp_c2 = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
            r2 = _rp_rank_profile_naive(mrp_r2, mrp_c2, A);
        }

        if (r2 != rank)
            res = 10;
        for (k = 0; k < rank && res == 0; k++)
            if (mrp_r2[k] != P[k] || mrp_c2[k] != Q[k])
                res = 11;

        if (mrp_r == NULL)
        {
            flint_free(mrp_r2);
            flint_free(mrp_c2);
        }
    }

    return res;
}

/*
    The basis of the right kernel of A with identity at the columns
    without pivot, from the reduced row echelon form of A computed by
    nmod_mat_rref (independently of nmod_mat_ple). X is initialised by
    this function, n x (n - rank); returns the rank.
*/
static slong
_rp_kernel_rref(nmod_mat_t X, const nmod_mat_t A)
{
    slong n = A->c, rank, i, j, k;
    slong * piv = flint_malloc((n + 1) * sizeof(slong));
    slong * nonpiv = flint_malloc((n + 1) * sizeof(slong));
    nmod_mat_t R;

    nmod_mat_init_set(R, A);
    rank = nmod_mat_rref(R);

    for (i = j = k = 0; j < n; j++)
    {
        if (i < rank && nmod_mat_entry(R, i, j) != 0)
            piv[i++] = j;
        else
            nonpiv[k++] = j;
    }

    nmod_mat_init(X, n, n - rank, A->mod.n);
    for (k = 0; k < n - rank; k++)
    {
        nmod_mat_entry(X, nonpiv[k], k) = 1;
        for (i = 0; i < rank; i++)
            nmod_mat_entry(X, piv[i], k) = nmod_neg(nmod_mat_entry(R, i, nonpiv[k]), A->mod);
    }

    nmod_mat_clear(R);
    flint_free(piv);
    flint_free(nonpiv);
    return rank;
}

/* random prime modulus, with the primes at the boundaries of the
   reductions of the delayed algorithms over-represented */
static ulong
_rp_rand_prime(flint_rand_t state)
{
    switch (n_randint(state, 12))
    {
        case 0:
            return 2;
        case 1:
            return 3;
        case 2:
            return n_nextprime(UWORD(1) << (FLINT_BITS - 1), 1);
        case 3:
#if FLINT_BITS == 64
            return UWORD(18446744073709551557);  /* largest prime < 2^64 */
#else
            return UWORD(4294967291);  /* largest prime < 2^32 */
#endif
        case 4:
#if FLINT_BITS == 64
            return UWORD(4294967291);  /* largest prime < 2^32 */
#else
            return UWORD(65521);  /* largest prime < 2^16 */
#endif
        case 5:
#if FLINT_BITS == 64
            return UWORD(9223372036854775783);  /* largest prime < 2^63 */
#else
            return UWORD(2147483647);  /* largest prime < 2^31 */
#endif
        case 6:
            return n_nextprime(UWORD(1) << (FLINT_BITS - 2), 1);
        default:
            return n_randtest_prime(state, 0);
    }
}

#endif
