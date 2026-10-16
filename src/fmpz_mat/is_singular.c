/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "perm.h"
#include "fmpz_mat.h"

/*
    Try to prove that the m x n matrix A has rank less than n by exhibiting
    a nonzero integer kernel vector, starting from the rank structure of A
    modulo a prime: rank, a row permutation P whose first rank entries index
    rows that are independent modulo p, and the columns pivs (pivot columns
    first, the remaining columns from position rank on), as produced by
    nmod_mat_lu_with_pivots. Requires rank < n.

    The pivot submatrix A[P[0..rank), pivs[0..rank)] is nonsingular modulo p,
    hence nonsingular over Q, and solving it against the first non-pivot
    column gives a candidate vector v (den at the non-pivot column, the
    negated solution at the pivot columns) with A[P[i]] v = 0 for i < rank by
    construction. If the remaining rows also annihilate v -- an exact check
    -- then v is a nonzero kernel vector and rank(A) < n unconditionally,
    whatever the prime. If the true rank exceeds the rank modulo p the check
    fails; it cannot succeed by accident.

    Compared with certifying the full rank as in fmpz_mat_rref_mul, this
    needs one right-hand side instead of n - rank of them, which matters
    when the corank is larger than one; for solving and for determinants,
    "singular" is all that is needed.

    Returns 1 if singularity was certified, 0 if the certificate failed
    (which for random primes has negligible probability when A is really
    rank deficient).
*/
int
_fmpz_mat_certify_singular_lu_mod_p(const fmpz_mat_t A,
    slong rank, const slong * P, const slong * pivs, ulong p)
{
    slong m = A->r, n = A->c, i, k, jfree;
    fmpz_mat_t As, rhs, z;
    fmpz_t den, acc, t;
    int ok;

    if (rank >= n)
        return 0;

    jfree = pivs[rank];

    if (rank == 0)
    {
        /* v = e_jfree: certifies iff column jfree of A is exactly zero */
        for (i = 0; i < m; i++)
            if (!fmpz_is_zero(fmpz_mat_entry(A, i, jfree)))
                return 0;
        return 1;
    }

    fmpz_mat_init(As, rank, rank);
    fmpz_mat_init(rhs, rank, 1);
    fmpz_mat_init(z, rank, 1);
    fmpz_init(den);
    fmpz_init(acc);
    fmpz_init(t);

    for (i = 0; i < rank; i++)
    {
        for (k = 0; k < rank; k++)
            fmpz_set(fmpz_mat_entry(As, i, k), fmpz_mat_entry(A, P[i], pivs[k]));
        fmpz_neg(fmpz_mat_entry(rhs, i, 0), fmpz_mat_entry(A, P[i], jfree));
    }

    /*
        As is the pivot submatrix, so it is nonsingular modulo p and the
        same prime can be handed to the p-adic solver together with the
        factorisation of As modulo p. This matters beyond saving a prime
        search: going through fmpz_mat_solve would let the solver's own
        prime search certify singularity, calling this function again on a
        smaller submatrix. Such a nested certificate always fails (As is
        nonsingular), but each failure would spawn another solve, so the
        work would grow like (bad primes per level)^depth -- readily
        reproducible with the small-primes testing hook. Forwarding p
        removes the recursion entirely.

        Small systems call the fraction-free algorithms directly rather
        than fmpz_mat_solve, whose dispatch could send them to a modular
        solver -- and hence back here -- if its cutoffs ever change.
    */
    if (rank <= 3)
        ok = fmpz_mat_solve_cramer(z, den, As, rhs);
    else if (rank <= 15)
        ok = fmpz_mat_solve_fflu(z, den, As, rhs);
    else
    {
        nmod_mat_t LU;
        slong * Ps = _perm_init(rank);

        nmod_mat_init(LU, rank, rank, p);
        fmpz_mat_get_nmod_mat(LU, As);
        ok = (nmod_mat_lu(Ps, LU, 0) == rank);
        if (ok)
            ok = _fmpz_mat_solve_dixon_den_given_lu(z, den, As, rhs, LU, Ps, p);

        nmod_mat_clear(LU);
        _perm_clear(Ps);
    }

    /* v has den at column jfree and z at the pivot columns; the rows
       P[0..rank) annihilate v by construction of the exact solution, so
       only the remaining rows need checking */
    for (i = rank; i < m && ok; i++)
    {
        fmpz_mul(acc, den, fmpz_mat_entry(A, P[i], jfree));
        for (k = 0; k < rank; k++)
        {
            fmpz_mul(t, fmpz_mat_entry(z, k, 0), fmpz_mat_entry(A, P[i], pivs[k]));
            fmpz_add(acc, acc, t);
        }
        ok = fmpz_is_zero(acc);
    }

    fmpz_mat_clear(As);
    fmpz_mat_clear(rhs);
    fmpz_mat_clear(z);
    fmpz_clear(den);
    fmpz_clear(acc);
    fmpz_clear(t);

    return ok;
}


static ulong
_det_3x3(ulong a, ulong b, ulong c, ulong d, ulong e, ulong f,
            ulong g, ulong h, ulong i)
{
    return a * (e * i - f * h) + b * (f * g - d * i) + c * (d * h - e * g);
}

static ulong
_det_4x4(const ulong * mm)
{
    ulong s, t, u, v;
#define m(ii,jj) mm[(ii) * 4 + (jj)]
    s = _det_3x3(m(1,1), m(1,2), m(1,3), m(2,1), m(2,2), m(2,3), m(3,1), m(3,2), m(3,3));
    t = _det_3x3(m(1,0), m(1,2), m(1,3), m(2,0), m(2,2), m(2,3), m(3,0), m(3,2), m(3,3));
    u = _det_3x3(m(1,0), m(1,1), m(1,3), m(2,0), m(2,1), m(2,3), m(3,0), m(3,1), m(3,3));
    v = _det_3x3(m(1,0), m(1,1), m(1,2), m(2,0), m(2,1), m(2,2), m(3,0), m(3,1), m(3,2));
    return m(0,0) * s - m(0,1) * t + m(0,2) * u - m(0,3) * v;
#undef m
}

/* Return x mod 2^FLINT_BITS */
static ulong _fmpz_get_ui(const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return (ulong) (*x);
    }
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*x);
        return (z->_mp_size > 0) ? z->_mp_d[0] : -z->_mp_d[0];
    }
}

int
fmpz_mat_is_singular(const fmpz_mat_t A)
{
    slong n = A->r, rank, rank_done = -1, bits_A;
    nmod_mat_t Am;
    slong * P, * pivs;
    fmpz_t det_bound, tested;
    ulong p;
    int result = 0;

    if (A->r != A->c)
        flint_throw(FLINT_ERROR, "non-square matrix in fmpz_mat_is_singular\n");

    if (n == 0)
        return 0;

    if (n == 1)
        return fmpz_is_zero(fmpz_mat_entry(A, 0, 0));

    if (n == 2)
    {
        ulong a = _fmpz_get_ui(fmpz_mat_entry(A, 0, 0));
        ulong b = _fmpz_get_ui(fmpz_mat_entry(A, 0, 1));
        ulong c = _fmpz_get_ui(fmpz_mat_entry(A, 1, 0));
        ulong d = _fmpz_get_ui(fmpz_mat_entry(A, 1, 1));

        /* det(A) != 0 mod 2^FLINT_BITS */
        if (a * d != b * c)
            return 0;

        if (!COEFF_IS_MPZ(*fmpz_mat_entry(A, 0, 0)) &&
            !COEFF_IS_MPZ(*fmpz_mat_entry(A, 0, 1)) &&  
            !COEFF_IS_MPZ(*fmpz_mat_entry(A, 1, 0)) &&  
            !COEFF_IS_MPZ(*fmpz_mat_entry(A, 1, 1)))
        {
            ulong u1, u0, v1, v0;
            smul_ppmm(u1, u0, a, d);
            smul_ppmm(v1, v0, b, c);
            return (u1 == v1);
        }

        fmpz_t u, v;
        fmpz_init(u);
        fmpz_init(v);
        fmpz_mul(u, fmpz_mat_entry(A, 0, 0), fmpz_mat_entry(A, 1, 1));
        fmpz_mul(v, fmpz_mat_entry(A, 0, 1), fmpz_mat_entry(A, 1, 0));
        result = fmpz_equal(u, v);
        fmpz_clear(u);
        fmpz_clear(v);
        return result;
    }

    if (n <= 4)
    {
        slong i, j;
        ulong m[16], d;
        fmpz_t det;

        /* det(A) mod 2^FLINT_BITS certifies det(A) != 0 at the cost of a
           few word operations, which settles the nonsingular case (for a
           singular matrix, and with probability 2^-FLINT_BITS otherwise,
           it is inconclusive and the determinant is computed exactly) */
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                m[i * n + j] = _fmpz_get_ui(fmpz_mat_entry(A, i, j));

        if (n == 3)
            d = _det_3x3(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
        else
            d = _det_4x4(m);

        if (d != 0)
            return 0;

        fmpz_init(det);
        fmpz_mat_det(det, A);
        result = fmpz_is_zero(det);
        fmpz_clear(det);
        return result;
    }

    /* the determinant bound is only needed to bound the number of bad
       primes, an event of negligible probability, so it is computed
       lazily */
    fmpz_init(det_bound);
    fmpz_init(tested);
    fmpz_one(tested);

    bits_A = FLINT_ABS(fmpz_mat_max_bits(A));

    P = flint_malloc(sizeof(slong) * n);
    pivs = flint_malloc(sizeof(slong) * n);
    nmod_mat_init(Am, n, n, 2);

    /* a small prime keeps the modular rank computation fast; the
       probability that it underestimates the rank is small, and a
       failure is caught by the certificate below (fmpz_mat_rref_mul
       makes the same choice) */
    p = flint_fmpz_mat_force_small_primes ? UWORD(1) : (UWORD(1) << 16);

    while (1)
    {
        p = n_nextprime(p, 1);
        nmod_mat_set_mod(Am, p);
        fmpz_mat_get_nmod_mat(Am, A);
        rank = nmod_mat_lu_with_pivots(P, pivs, Am);

        /* full rank modulo p certifies nonsingularity outright */
        if (rank == n)
            break;

        /* Rank deficient modulo p, so the matrix is probably singular.
           For a small matrix with a small determinant bound, computing
           the determinant outright is cheaper than the kernel-vector
           certificate, which involves an exact solve of comparable size
           (measured crossover; the certificate wins from about n = 16,
           and by orders of magnitude for large n). */
        if (n <= 16 && (double) n * (bits_A + 4) <= 12000.0)
        {
            fmpz_t det;
            fmpz_init(det);
            fmpz_mat_det(det, A);
            result = fmpz_is_zero(det);
            fmpz_clear(det);
            break;
        }

        /* The certificate succeeds exactly when the rank modulo p is the
           true rank, so a failure at rank r proves the true rank exceeds
           r and no prime of rank at most r can succeed either: attempt it
           only when the rank increases. */
        if (rank > rank_done)
        {
            rank_done = rank;

            if (_fmpz_mat_certify_singular_lu_mod_p(A, rank, P, pivs, p))
            {
                result = 1;
                break;
            }
        }

        /* bad prime (rank modulo p below the true rank): such a prime
           divides some maximal minor, so once the product of the tested
           primes exceeds the determinant bound, and every tested prime
           saw a rank drop, the matrix really is singular */
        if (fmpz_is_zero(det_bound))
        {
            fmpz_mat_det_bound(det_bound, A);
            if (fmpz_is_zero(det_bound))
            {
                /* a zero row or column */
                result = 1;
                break;
            }
        }

        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
        {
            result = 1;
            break;
        }
    }

    nmod_mat_clear(Am);
    flint_free(P);
    flint_free(pivs);
    fmpz_clear(det_bound);
    fmpz_clear(tested);

    return result;
}
