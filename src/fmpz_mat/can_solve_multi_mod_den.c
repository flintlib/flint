/*
    Copyright (C) 2020 William Hart
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "perm.h"
#include "thread_support.h"


/*
    Solve A X = B, A being m x n of rank r (determined mod p). We pick the
    r pivot rows and pivot columns of A, solve the nonsingular r x r
    system for the pivot-column entries of X (other entries zero) using
    the given square solver (Dixon lifting or multimodular), and then verify the remaining m - r
    equations exactly. If the verification fails, either the rank mod p
    was too small or the system is inconsistent; we try more primes until
    the product of primes exceeds the determinant bound (which certifies
    that the rank mod p was correct at least once).
*/
typedef int (*square_solver_fn)(fmpz_mat_t, fmpz_t, const fmpz_mat_t, const fmpz_mat_t);

typedef struct
{
    nmod_mat_struct * LUb;
    slong ** permb;
    slong * rankb;
    slong m;
} cs_lu_args;

static void
cs_lu_worker(slong k, void * arg)
{
    cs_lu_args * a = (cs_lu_args *) arg;
    slong i;
    for (i = 0; i < a->m; i++)
        a->permb[k][i] = i;
    a->rankb[k] = nmod_mat_lu(a->permb[k], a->LUb + k, 0);
}

static int
_fmpz_mat_can_solve_generic(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B, square_solver_fn solve)
{
    ulong p;
    fmpz_t tested;
    nmod_mat_t LU;
    int result = 0, success;
    slong * perm = NULL, * pivots;
    slong i, j, k, col, rank, rank_done = 0, m, n, cols;
    int certified_inconsistent = 0;
    /* primes are processed in growing batches: A is reduced modulo the
       whole batch at once and the rank computations run in parallel */
#define CS_BATCH 32
    nmod_mat_t LUb[CS_BATCH];
    slong * permb[CS_BATCH];
    slong rankb[CS_BATCH];
    ulong pb[CS_BATCH];
    slong nbatch = 0, ibatch = 0, batch_size = 1, nalloc_batch = 0, batch_limit;
    fmpz_mat_t Arank, Brank, Zrank, Anp, Bnp, T;
    fmpz_t det_bound;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_can_solve). Incompatible matrix dimensions.\n");
    }

    m = A->r;
    n = A->c;
    cols = B->c;

    if (m == 0 || cols == 0)
    {
        fmpz_mat_zero(X);
        fmpz_one(den);
        return 1;
    }

    if (n == 0)
    {
        fmpz_mat_zero(X);
        fmpz_one(den);
        return fmpz_mat_is_zero(B);
    }

    /* X = 0 solves A X = 0 whatever the rank of A */
    if (fmpz_mat_is_zero(B))
    {
        fmpz_mat_zero(X);
        fmpz_one(den);
        return 1;
    }

    p = (flint_fmpz_mat_force_small_primes ? UWORD(1) : (UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS));

    fmpz_init(det_bound);
    fmpz_init(tested);
    fmpz_one(tested);

    nmod_mat_init(LU, m, n, p);

    pivots = flint_malloc(sizeof(slong) * n);

    /* Bound on all r x r minors of A: a prime is bad only if it divides
       some r x r minor, so once the product of the tested primes exceeds
       this bound at least one tested prime gave the true rank. */
    fmpz_mat_det_bound_submatrix(det_bound, A);

    /* batched reduction pays off for multi-limb entries; for word-size
       entries (one or two limbs) a plain per-prime reduction is cheaper */
    batch_limit = (FLINT_ABS(fmpz_mat_max_bits(A)) >= 64) ? CS_BATCH : 1;

    while (1)
    {
        if (ibatch == nbatch)
        {
            cs_lu_args largs;
            for (k = nalloc_batch; k < batch_size; k++)
            {
                nmod_mat_init(LUb[k], m, n, 2);
                permb[k] = _perm_init(m);
            }
            nalloc_batch = FLINT_MAX(nalloc_batch, batch_size);
            for (k = 0; k < batch_size; k++)
            {
                p = n_nextprime(p, 0);
                pb[k] = p;
                nmod_mat_set_mod(LUb[k], p);
            }
            if (batch_size == 1)
                fmpz_mat_get_nmod_mat(LUb[0], A);
            else
                fmpz_mat_multi_mod_ui(LUb, batch_size, A);
            largs.LUb = (nmod_mat_struct *) LUb;
            largs.permb = permb;
            largs.rankb = rankb;
            largs.m = m;
            flint_parallel_do(cs_lu_worker, &largs, batch_size,
                FLINT_MIN(FLINT_MAX(1, (slong) ((double) batch_size * m * n * FLINT_MIN(m, n) / 3 / FMPZ_MAT_SOLVE_MIN_WORK_PER_THREAD)), 1024),
                FLINT_PARALLEL_UNIFORM);
            nbatch = batch_size;
            ibatch = 0;
            batch_size = FLINT_MIN(batch_limit, 2 * batch_size);
        }
        k = ibatch++;
        p = pb[k];
        rank = rankb[k];
        perm = permb[k];
        nmod_mat_swap(LU, LUb[k]);
        nmod_mat_set_mod(LUb[k], 2);

        /* B is nonzero, so a rank 0 matrix cannot solve the system; but
           the rank mod p may be too small, so certification (below) is
           still needed */
        if (rank == 0)
            goto bad_prime;

        /* The exact candidate solution and its verification depend only
           on the rank: any prime whose rank equals the true rank is a
           good prime, so all primes attaining the maximal rank seen give
           the same (correct, once one of them is good) verdict. Exact
           work is therefore redone only when the rank increases; primes
           with rank <= rank_done need just this nmod LU, which makes
           certifying an inconsistent system cheap. */
        if (rank <= rank_done)
            goto bad_prime;
        rank_done = rank;

        /* square and nonsingular mod p: the system has a unique candidate
           solution and every equation is a pivot equation */
        if (rank == m && rank == n)
        {
            result = solve(X, den, A, B);
            break;
        }

        col = 0;
        for (i = 0; i < rank; i++)
        {
            while (nmod_mat_entry(LU, i, col) == 0)
                col++;
            pivots[i] = col;
            col++;
        }

        fmpz_mat_init(Arank, rank, rank);
        fmpz_mat_init(Brank, rank, cols);
        fmpz_mat_init(Zrank, rank, cols);
        fmpz_mat_init(Anp, m - rank, rank);
        fmpz_mat_init(Bnp, m - rank, cols);
        fmpz_mat_init(T, m - rank, cols);

        for (i = 0; i < rank; i++)
        {
            for (k = 0; k < rank; k++)
                fmpz_set(fmpz_mat_entry(Arank, i, k), fmpz_mat_entry(A, perm[i], pivots[k]));
            for (j = 0; j < cols; j++)
                fmpz_set(fmpz_mat_entry(Brank, i, j), fmpz_mat_entry(B, perm[i], j));
        }
        for (i = rank; i < m; i++)
        {
            for (k = 0; k < rank; k++)
                fmpz_set(fmpz_mat_entry(Anp, i - rank, k), fmpz_mat_entry(A, perm[i], pivots[k]));
            for (j = 0; j < cols; j++)
                fmpz_set(fmpz_mat_entry(Bnp, i - rank, j), fmpz_mat_entry(B, perm[i], j));
        }

        success = solve(Zrank, den, Arank, Brank);

        if (success)
        {
            /* the r pivot equations hold by construction; check the rest */
            fmpz_mat_mul(T, Anp, Zrank);
            fmpz_mat_scalar_mul_fmpz(Bnp, Bnp, den);
            result = fmpz_mat_equal(T, Bnp);

            if (result)
            {
                fmpz_mat_zero(X);
                for (i = 0; i < rank; i++)
                    for (j = 0; j < cols; j++)
                        fmpz_swap(fmpz_mat_entry(X, pivots[i], j), fmpz_mat_entry(Zrank, i, j));
            }
            else
            {
                /* The check failed at some non-pivot row i0. Try to prove
                   inconsistency outright with a certificate of
                   inconsistency in the sense of the Fredholm alternative
                   (a vector y with y^t A = 0 and y^t B != 0; see
                   [GieLobSau1998] for its use as a computational
                   certificate and [MulSto2004] for certified dense
                   solving over Q): solve the transpose system
                   (A_pivot)^t z = (row i0 of A restricted to the pivot
                   columns)^t exactly and let
                   y = wden e_{i0} - sum_k z_k e_{perm[k]}. Then y^t A
                   vanishes on the pivot columns by construction, and
                   y^t B != 0 holds automatically because row i0 failed
                   the check above. If y^t A also vanishes on the
                   remaining columns -- an exact, self-certifying
                   computation -- then no X can satisfy A X = den' B, and
                   the system is proven inconsistent with no assumption
                   on the prime. Otherwise the rank mod p was too small
                   and we fall through to the bad-prime loop, which is an
                   event of negligible probability for random primes. */
                slong i0 = -1;

                for (i = rank; i < m && i0 < 0; i++)
                    for (j = 0; j < cols; j++)
                        if (!fmpz_equal(fmpz_mat_entry(T, i - rank, j),
                                        fmpz_mat_entry(Bnp, i - rank, j)))
                        {
                            i0 = i;
                            break;
                        }

                if (i0 >= 0)
                {
                    fmpz_mat_t At, bt, z;
                    fmpz_t wden, acc, t2;
                    int witness;

                    fmpz_mat_init(At, rank, rank);
                    fmpz_mat_init(bt, rank, 1);
                    fmpz_mat_init(z, rank, 1);
                    fmpz_init(wden);
                    fmpz_init(acc);
                    fmpz_init(t2);

                    fmpz_mat_transpose(At, Arank);
                    for (k = 0; k < rank; k++)
                        fmpz_set(fmpz_mat_entry(bt, k, 0), fmpz_mat_entry(A, perm[i0], pivots[k]));

                    witness = solve(z, wden, At, bt);

                    /* y^t A on all columns: wden A[i0] - z^t A[perm rows] */
                    for (j = 0; j < n && witness; j++)
                    {
                        fmpz_mul(acc, wden, fmpz_mat_entry(A, perm[i0], j));
                        for (k = 0; k < rank; k++)
                        {
                            fmpz_mul(t2, fmpz_mat_entry(z, k, 0), fmpz_mat_entry(A, perm[k], j));
                            fmpz_sub(acc, acc, t2);
                        }
                        witness = fmpz_is_zero(acc);
                    }

                    /* y^t B != 0 (automatic; checked for safety) */
                    if (witness)
                    {
                        int nz = 0;
                        for (j = 0; j < cols && !nz; j++)
                        {
                            fmpz_mul(acc, wden, fmpz_mat_entry(B, perm[i0], j));
                            for (k = 0; k < rank; k++)
                            {
                                fmpz_mul(t2, fmpz_mat_entry(z, k, 0), fmpz_mat_entry(B, perm[k], j));
                                fmpz_sub(acc, acc, t2);
                            }
                            nz = !fmpz_is_zero(acc);
                        }
                        witness = nz;
                    }

                    fmpz_mat_clear(At);
                    fmpz_mat_clear(bt);
                    fmpz_mat_clear(z);
                    fmpz_clear(wden);
                    fmpz_clear(acc);
                    fmpz_clear(t2);

                    certified_inconsistent = witness;
                }
            }
        }

        fmpz_mat_clear(Arank);
        fmpz_mat_clear(Brank);
        fmpz_mat_clear(Zrank);
        fmpz_mat_clear(Anp);
        fmpz_mat_clear(Bnp);
        fmpz_mat_clear(T);

        if (result || certified_inconsistent)
            break;
bad_prime: ;

        /* A prime p is "bad" only if the rank mod p is less than the
           true rank, which requires p to divide some r x r minor; once
           the product of tested primes exceeds the bound on the minors,
           at least one prime gave the true rank, so the system really is
           inconsistent. */
        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
            break;
    }

    if (!result)
        fmpz_one(den);

    fmpz_clear(det_bound);
    for (k = 0; k < nalloc_batch; k++)
    {
        nmod_mat_clear(LUb[k]);
        _perm_clear(permb[k]);
    }
    nmod_mat_clear(LU);
    fmpz_clear(tested);
    flint_free(pivots);

    return result;
}

#undef CS_BATCH

int
fmpz_mat_can_solve_dixon_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_can_solve_generic(X, den, A, B, fmpz_mat_solve_dixon_den);
}

int
fmpz_mat_can_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_can_solve_generic(X, den, A, B, fmpz_mat_solve_multi_mod_den);
}

/* the square subsystem is solved with fmpz_mat_solve, which picks the
   solver by shape */
int
_fmpz_mat_can_solve_auto_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    return _fmpz_mat_can_solve_generic(X, den, A, B, fmpz_mat_solve);
}
