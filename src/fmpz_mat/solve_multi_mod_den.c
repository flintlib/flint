/*
    Copyright (C) 2011, 2026 Fredrik Johansson
    Copyright (C) 2020 William Hart
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "perm.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "thread_support.h"

/*
    Find a prime p such that A is nonsingular mod p, and solve A X = B
    mod p. Returns p, or 0 if A is certified singular. On success the
    solution is left in Xmod (with modulus p).
*/
static ulong
_find_good_prime_and_solve(nmod_mat_t Xmod, nmod_mat_t Amod, nmod_mat_t Bmod,
        const fmpz_mat_t A, const fmpz_mat_t B, const fmpz_t det_bound)
{
    ulong p;
    slong i, n, rank, rank_done = -1, * pivs, * P;
    fmpz_t tested;
    nmod_mat_t PB;

    n = A->r;
    pivs = (slong *) flint_malloc(n * sizeof(slong));
    P = _perm_init(n);
    fmpz_init(tested);
    fmpz_one(tested);

    p = (flint_fmpz_mat_force_small_primes ? UWORD(1) : (UWORD(1) << NMOD_MAT_OPTIMAL_MODULUS_BITS));

    while (1)
    {
        p = n_nextprime(p, 0);
        nmod_mat_set_mod(Xmod, p);
        nmod_mat_set_mod(Amod, p);
        nmod_mat_set_mod(Bmod, p);
        fmpz_mat_get_nmod_mat(Amod, A);
        fmpz_mat_get_nmod_mat(Bmod, B);
        rank = nmod_mat_lu_with_pivots(P, pivs, Amod);

        if (rank == n)
        {
            nmod_mat_init(PB, B->r, B->c, p);
            for (i = 0; i < n; i++)
                _nmod_vec_set(nmod_mat_entry_ptr(PB, i, 0),
                              nmod_mat_entry_ptr(Bmod, P[i], 0), Bmod->c);
            nmod_mat_solve_tril(Xmod, Amod, PB, 1);
            nmod_mat_solve_triu(Xmod, Amod, Xmod, 0);
            nmod_mat_clear(PB);
            break;
        }

        fmpz_mul_ui(tested, tested, p);
        if (fmpz_cmp(tested, det_bound) > 0)
        {
            p = 0;
            break;
        }

        /* The kernel-vector certificate succeeds exactly when the rank
           modulo p is the true rank, so a failure at rank r rules out
           every prime of rank at most r; attempt it only when the rank
           increases. */
        if (rank > rank_done)
        {
            rank_done = rank;

            if (_fmpz_mat_certify_singular_lu_mod_p(A, rank, P, pivs, p))
            {
                p = 0;
                break;
            }
        }
    }

    fmpz_clear(tested);
    flint_free(pivs);
    _perm_clear(P);
    return p;
}

typedef struct
{
    nmod_mat_struct * Abatch;
    nmod_mat_struct * Bbatch;
    nmod_mat_struct * Xbatch;
    int * ok;
} mm_solve_args;

static void
mm_solve_worker(slong k, void * arg)
{
    mm_solve_args * a = (mm_solve_args *) arg;
    a->ok[k] = nmod_mat_solve(a->Xbatch + k, a->Abatch + k, a->Bbatch + k);
}

/*
    Multimodular solving of A X = B with A nonsingular. Residues of the
    solution are accumulated for a growing set of primes and combined with
    fmpz_mat_multi_CRT_ui only when a reconstruction attempt is made, so the
    total CRT cost is a small multiple of one full CRT rather than
    quadratic in the number of primes. Termination is output-sensitive:
    a random projection of the solution is reconstructed cheaply for each
    new prime, and when it stabilises the full solution is reconstructed
    and certified (see _fmpz_mat_solve_reconstruct_attempt).
*/
static int
_fmpz_mat_solve_multi_mod_den(fmpz_mat_t Z, fmpz_t den,
        const fmpz_mat_t A, const fmpz_mat_t B,
        nmod_mat_t Xmod, ulong p, const fmpz_t N, const fmpz_t D)
{
    int found = 1;
    fmpz_t bound, pprod, Amax, Bmax, s, sprod, snum, sden, snum_prev, sden_prev, t;
    /* attempts reconstruct into Zt, swapped into Z only on success, so
       that Z may alias A or B (as fmpz_mat_inv does) */
    fmpz_mat_t x, Zt;
    nmod_mat_t * res;
    slong alloc, i, j, k, n, cols, next_full, last_probe, probe_gap;
    /* A and B are reduced modulo batches of primes at once (much cheaper
       than one prime at a time for large entries) */
#define MM_BATCH 32
    nmod_mat_t Abatch[MM_BATCH], Bbatch[MM_BATCH], Xbatch[MM_BATCH];
    int okbatch[MM_BATCH];
    ulong pbatch[MM_BATCH];
    slong nbatch = 0, ibatch = 0, batch_size, nalloc_batch = 0;
    ulong * u;
    ulong sp;
    slong den_margin = 16;
    int probe_ok = 0, probe_stable = 0, x_valid = 0;
    flint_rand_t state;

    n = A->r;
    cols = B->c;

    fmpz_init(bound);
    fmpz_init(pprod);
    fmpz_init(Amax);
    fmpz_init(Bmax);
    fmpz_init(s);
    fmpz_init(sprod);
    fmpz_init(snum);
    fmpz_init(sden);
    fmpz_init(snum_prev);
    fmpz_init(sden_prev);
    fmpz_init(t);
    fmpz_mat_init(x, n, cols);
    fmpz_mat_init(Zt, n, cols);

    /* modulus needed for guaranteed reconstruction: m > 2 N D */
    fmpz_mul(bound, N, D);
    fmpz_mul_2exp(bound, bound, 1);

    /* track pointers to the largest entries; only two multiprecision
       operations at the end */
    {
        const fmpz * amax = fmpz_mat_entry(A, 0, 0);
        const fmpz * bmax = fmpz_mat_entry(B, 0, 0);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                if (fmpz_cmpabs(amax, fmpz_mat_entry(A, i, j)) < 0)
                    amax = fmpz_mat_entry(A, i, j);
        for (i = 0; i < n; i++)
            for (j = 0; j < cols; j++)
                if (fmpz_cmpabs(bmax, fmpz_mat_entry(B, i, j)) < 0)
                    bmax = fmpz_mat_entry(B, i, j);
        fmpz_abs(Amax, amax);
        fmpz_mul_ui(Amax, Amax, n);
        fmpz_abs(Bmax, bmax);
    }

    /* random projection u . x[:,0], maintained as an fmpz by incremental
       CRT (a single number; its cost is negligible next to the solves) */
    flint_rand_init(state);
    u = (ulong *) flint_malloc(sizeof(ulong) * n * cols);
    for (i = 0; i < n * cols; i++)
        u[i] = n_randint(state, UWORD(1) << 4) + 1;

    alloc = 16;
    res = (nmod_mat_t *) flint_malloc(sizeof(nmod_mat_t) * alloc);
    nmod_mat_init(res[0], n, cols, p);
    nmod_mat_set(res[0], Xmod);
    fmpz_set_ui(pprod, p);
    sp = 0;
    /* the weights are < 2^4, so they are already reduced except for the
       tiny primes used by the small-primes testing hook */
    for (k = 0; k < n; k++)
        for (j = 0; j < cols; j++)
            sp = nmod_add(sp, nmod_mul(nmod_mat_entry(res[0], k, j),
                p > 16 ? u[k * cols + j] : u[k * cols + j] % p, res[0]->mod), res[0]->mod);
    fmpz_set_ui(s, sp);
    fmpz_set_ui(sprod, p);

    i = 1; /* number of primes */
    /* the batch size starts at 1 and doubles with each batch, so that at
       most about half of the primes reduced and solved are wasted when the
       solution is found early, while long lifts soon use full batches */
    batch_size = 1;
    next_full = 1;
    last_probe = 0;
    probe_gap = 1;

    while (1)
    {
        if (fmpz_cmp(pprod, bound) > 0)
        {
            fmpz_mat_multi_CRT_ui(x, res, i, 1);
            fmpz_one(t);
            /* The reconstruction is unique for m > 2 N D and the true
               solution satisfies the bounds, so this cannot fail; report
               failure rather than a wrong answer if it ever does. */
            if (!_fmpz_mat_reconstruct_matwise(Zt, den, x, pprod, N, D, t))
            {
                found = 0;
                goto done;
            }
            fmpz_mat_swap(Z, Zt);
            goto done;
        }

        if (i - last_probe >= probe_gap)
        {
            last_probe = i;
            probe_ok = _fmpz_mat_solve_probe_reconstruct(snum, sden, s, pprod, sden_prev, probe_ok);
            probe_gap = FLINT_MAX(i / 16, 1);
            /* a successful probe carries a 2^-30 confidence margin (see
               _fmpz_mat_solve_probe_reconstruct); no need to wait for two
               agreeing probes */
            if (probe_ok)
                probe_stable = 1;
            else
            {
                probe_stable = 0;
                next_full = i;
            }
            fmpz_swap(snum, snum_prev);
            fmpz_swap(sden, sden_prev);
        }

        /* while the modulus is small an attempt is cheap and the probe (a
           weighted sum needing a larger modulus than the entries, and
           often spuriously successful) is not informative: attempt blindly
           with balanced bounds and rely on the certificate / verification */
        if (((probe_ok && probe_stable) || i <= 8) && i >= next_full
            && (i <= 8 || _fmpz_mat_solve_attempt_worthwhile(pprod, snum_prev, sden_prev, Amax, Bmax, n, cols,
                    NMOD_MAT_OPTIMAL_MODULUS_BITS - 1,
                    (2.0 * (double) (n / 3 + cols) * 64.0 / ((double) cols * (1 + fmpz_bits(Amax))) + 1.0 / cols))))
        {
            int blind = !probe_ok || i <= 8;
            if (!x_valid)
            {
                fmpz_mat_multi_CRT_ui(x, res, i, 1);
                x_valid = 1;
            }
            /* a step (one more prime) costs an nmod LU + solve, i.e. about
               (n/3 + cols) n^2 modular operations (each ~2 word
               operations of the unit), plus the reduction of A
               modulo the prime (~ n^2 bits(A)/64), versus n^2 cols
               bits(A)/64 for multiplying A by a word-size matrix */
            if (blind)
                fmpz_one(t);
            if (_fmpz_mat_solve_reconstruct_attempt(Zt, den, x, pprod, N, D,
                    blind ? t : sden_prev, Amax, Bmax, A, B, NMOD_MAT_OPTIMAL_MODULUS_BITS - 1,
                    (2.0 * (double) (n / 3 + cols) * 64.0 / ((double) cols * (1 + fmpz_bits(Amax))) + 1.0 / cols), blind ? 0 : den_margin))
            {
                fmpz_mat_swap(Z, Zt);
                goto done;
            }
            if (probe_stable && den_margin < (slong) fmpz_bits(D))
                den_margin *= 2;
            next_full = i + FLINT_MAX(1, i / 4);
        }

        /* next good prime, from a batch of primes reduced at once */
        while (1)
        {
            if (ibatch == nbatch)
            {
                mm_solve_args sargs;
                /* batch matrices are allocated as the batch grows */
                for (k = nalloc_batch; k < batch_size; k++)
                {
                    nmod_mat_init(Abatch[k], n, n, 2);
                    nmod_mat_init(Bbatch[k], n, cols, 2);
                    nmod_mat_init(Xbatch[k], n, cols, 2);
                }
                nalloc_batch = FLINT_MAX(nalloc_batch, batch_size);
                for (k = 0; k < batch_size; k++)
                {
                    p = n_nextprime(p, 1);
                    pbatch[k] = p;
                    nmod_mat_set_mod(Abatch[k], p);
                    nmod_mat_set_mod(Bbatch[k], p);
                    nmod_mat_set_mod(Xbatch[k], p);
                }
                fmpz_mat_multi_mod_ui(Abatch, batch_size, A);
                fmpz_mat_multi_mod_ui(Bbatch, batch_size, B);
                /* the modular solves of a batch are independent: run them
                   in parallel (with one thread the batch is as cheap as
                   solving prime by prime; the batch size grows
                   geometrically so that few solves are wasted when the
                   solution is found early) */
                sargs.Abatch = (nmod_mat_struct *) Abatch;
                sargs.Bbatch = (nmod_mat_struct *) Bbatch;
                sargs.Xbatch = (nmod_mat_struct *) Xbatch;
                sargs.ok = okbatch;
                /* thread_limit 1 makes flint_parallel_do a serial loop */
                flint_parallel_do(mm_solve_worker, &sargs, batch_size,
                    FLINT_MIN(FLINT_MAX(1, (slong) ((double) batch_size * n * n * (n / 3 + cols) / FMPZ_MAT_SOLVE_MIN_WORK_PER_THREAD)), 1024),
                    FLINT_PARALLEL_UNIFORM);
                nbatch = batch_size;
                ibatch = 0;
                /* grow the batch geometrically as the lifting goes on */
                batch_size = FLINT_MIN(MM_BATCH, 2 * batch_size);
            }
            k = ibatch++;
            if (okbatch[k])
            {
                nmod_mat_set_mod(Xmod, pbatch[k]);
                nmod_mat_swap(Xmod, Xbatch[k]);
                p = pbatch[k];
                break;
            }
        }

        if (i == alloc)
        {
            alloc *= 2;
            res = (nmod_mat_t *) flint_realloc(res, sizeof(nmod_mat_t) * alloc);
        }
        /* keep the residues (swap, no copy): the empty matrix given to
           res[i] ends up in Xmod, which is then given its proper size
           again (clearing it first, since an empty nmod_mat may still own
           an allocation) */
        nmod_mat_init(res[i], 0, 0, p);
        nmod_mat_swap(res[i], Xmod);
        nmod_mat_clear(Xmod);
        nmod_mat_init(Xmod, n, cols, p);
        x_valid = 0;

        sp = 0;
        for (k = 0; k < n; k++)
            for (j = 0; j < cols; j++)
                sp = nmod_add(sp, nmod_mul(nmod_mat_entry(res[i], k, j),
                    p > 16 ? u[k * cols + j] : u[k * cols + j] % p, res[i]->mod), res[i]->mod);
        fmpz_CRT_ui(s, s, sprod, sp, p, 0);
        fmpz_mul_ui(sprod, sprod, p);

        fmpz_mul_ui(pprod, pprod, p);
        i++;
    }

done:
    for (k = 0; k < nalloc_batch; k++)
    {
        nmod_mat_clear(Abatch[k]);
        nmod_mat_clear(Bbatch[k]);
        nmod_mat_clear(Xbatch[k]);
    }
    for (j = 0; j < i; j++)
        nmod_mat_clear(res[j]);
    flint_free(res);
    flint_free(u);
    flint_rand_clear(state);

    fmpz_clear(bound);
    fmpz_clear(pprod);
    fmpz_clear(Amax);
    fmpz_clear(Bmax);
    fmpz_clear(s);
    fmpz_clear(sprod);
    fmpz_clear(snum);
    fmpz_clear(sden);
    fmpz_clear(snum_prev);
    fmpz_clear(sden_prev);
    fmpz_clear(t);
    fmpz_mat_clear(x);
    fmpz_mat_clear(Zt);

    return found;
}

#undef MM_BATCH

int
fmpz_mat_solve_multi_mod_den(fmpz_mat_t X, fmpz_t den,
                        const fmpz_mat_t A, const fmpz_mat_t B)
{
    nmod_mat_t Xmod, Amod, Bmod;
    fmpz_t N, D;
    ulong p;

    if (!fmpz_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_solve_multi_mod_den). Non-square system matrix.\n");
    }

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
    {
        fmpz_one(den);
        return 1;
    }

    fmpz_init(N);
    fmpz_init(D);
    fmpz_mat_solve_bound(N, D, A, B);

    nmod_mat_init(Amod, A->r, A->c, 1);
    nmod_mat_init(Bmod, B->r, B->c, 1);
    nmod_mat_init(Xmod, B->r, B->c, 1);

    p = _find_good_prime_and_solve(Xmod, Amod, Bmod, A, B, D);
    if (p != 0)
    {
        if (!_fmpz_mat_solve_multi_mod_den(X, den, A, B, Xmod, p, N, D))
            p = 0;
    }

    nmod_mat_clear(Xmod);
    nmod_mat_clear(Bmod);
    nmod_mat_clear(Amod);
    fmpz_clear(N);
    fmpz_clear(D);

    return p != 0;
}
