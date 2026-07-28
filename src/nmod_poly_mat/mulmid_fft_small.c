/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/*
    Set each entry of C to the coefficient window [zl, zh) of the
    corresponding entry of the matrix product A*B, computed with a single
    fft_small plan: every input polynomial is transformed once, the k
    products of each output entry are accumulated pointwise in transform
    space, and one inverse transform plus chinese remaindering per output
    entry recovers the window.

    Work is distributed over the available threads at the task level:
    input transforms in one pass, then (accumulate, ifft, export) per
    output entry in a second pass. When there are fewer tasks than
    threads, the transform and export functions pick up the spare pool
    threads internally through their own flint_request_threads calls, so
    both small matrices of huge polynomials and large matrices of small
    polynomials parallelize.

    Returns 1 on success and 0 (leaving C undefined) if no plan exists,
    which requires the accumulation bound k*min(amax, bmax)*2^(2*modbits)
    to exceed the available primes and is unreachable for practical
    dimensions.
*/

typedef struct {
    fft_small_op_struct* ops;
    const nmod_poly_struct** tasks;
    ulong ntasks;
    nmod_t mod;
    const fft_small_plan_struct* P;
    ulong start;
    ulong stop;
} _transform_worker_struct;

static void _transform_worker_func(void* varg)
{
    _transform_worker_struct* W = (_transform_worker_struct*) varg;

    for (ulong t = W->start; t < W->stop; t++)
    {
        const nmod_poly_struct* p = W->tasks[t];

        if (p->length > 0)
            fft_small_fft_nmod(W->ops + t, p->coeffs, p->length,
                    n_round_up(p->length, BLK_SZ), W->mod, W->P);
    }
}

typedef struct {
    nmod_poly_mat_struct* C;
    const fft_small_op_struct* Aops;    /* entry (i, l) at index i*k + l */
    const fft_small_op_struct* Bops;    /* entry (l, j) at index l*bc + j */
    slong k;
    slong bc;
    ulong wlen;                         /* zh - zl of the original window */
    nmod_t mod;
    const fft_small_plan_struct* P;
    fft_small_op_struct Z[1];           /* per-worker accumulator */
    ulong* tmp;                         /* per-worker export buffer */
    ulong start_ij;
    ulong stop_ij;
} _entry_worker_struct;

static void _entry_worker_func(void* varg)
{
    _entry_worker_struct* W = (_entry_worker_struct*) varg;
    const fft_small_plan_struct* P = W->P;
    ulong elen = P->zh - P->zl;         /* exported length (window clamped
                                           to the maximal product length) */

    for (ulong idx = W->start_ij; idx < W->stop_ij; idx++)
    {
        slong i = idx / W->bc;
        slong j = idx % W->bc;
        nmod_poly_struct* entry = nmod_poly_mat_entry(W->C, i, j);
        ulong cnt = 0;

        for (slong l = 0; l < W->k; l++)
        {
            const fft_small_op_struct* Ail = W->Aops + i*W->k + l;
            const fft_small_op_struct* Blj = W->Bops + l*W->bc + j;

            /* zero factors were never transformed and contribute nothing */
            if (Ail->itrunc == 0 || Blj->itrunc == 0)
                continue;

            if (cnt == 0)
                fft_small_op_mul(W->Z, Ail, Blj, P);
            else
                fft_small_op_addmul(W->Z, Ail, Blj, P);
            cnt++;
        }

        nmod_poly_fit_length(entry, W->wlen);

        if (cnt == 0)
        {
            _nmod_vec_zero(entry->coeffs, W->wlen);
        }
        else
        {
            fft_small_ifft(W->Z, P);
            fft_small_export_nmod(W->tmp, W->Z, W->mod, P);
            flint_mpn_copyi(entry->coeffs, W->tmp, elen);
            if (W->wlen > elen)
                _nmod_vec_zero(entry->coeffs + elen, W->wlen - elen);
        }

        _nmod_poly_set_length(entry, W->wlen);
        _nmod_poly_normalise(entry);
    }
}

int nmod_poly_mat_mulmid_fft_small(nmod_poly_mat_t C,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                    slong zl, slong zh)
{
    slong ar = A->r;
    slong k = A->c;
    slong bc = B->c;
    slong i, j, l;
    ulong amax, bmax, zn, wlen, modbits, len_bound, hi;
    ulong nA, nB, nC, ntasks, t;
    nmod_t mod;
    fft_small_plan_t P;
    fft_small_op_struct* ops;
    const nmod_poly_struct** tasks;
    _transform_worker_struct* targs;
    _entry_worker_struct* eargs;
    thread_pool_handle* handles;
    slong nworkers;
    ulong nthreads;
    int success;

    FLINT_ASSERT(k == B->r);
    FLINT_ASSERT(C->r == ar && C->c == bc);
    FLINT_ASSERT(zl >= 0);
    FLINT_ASSERT(zh >= 0);

    if (ar == 0 || bc == 0)
        return 1;

    if (zl >= zh)
    {
        nmod_poly_mat_zero(C);
        return 1;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, ar, bc, nmod_poly_mat_modulus(A));
        success = nmod_poly_mat_mulmid_fft_small(T, A, B, zl, zh);
        if (success)
            nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return success;
    }

    nmod_init(&mod, nmod_poly_mat_modulus(A));
    wlen = (ulong) (zh - zl);

    amax = 0;
    for (i = 0; i < ar; i++)
        for (l = 0; l < k; l++)
            amax = FLINT_MAX(amax, (ulong) nmod_poly_length(nmod_poly_mat_entry(A, i, l)));
    bmax = 0;
    for (l = 0; l < k; l++)
        for (j = 0; j < bc; j++)
            bmax = FLINT_MAX(bmax, (ulong) nmod_poly_length(nmod_poly_mat_entry(B, l, j)));

    /* every product lies beyond the window, or there are no products */
    if (amax == 0 || bmax == 0 || (ulong) zl >= amax + bmax - 1)
    {
        nmod_poly_mat_zero(C);
        return 1;
    }

    zn = amax + bmax - 1;
    hi = FLINT_MIN((ulong) zh, zn);

    modbits = FLINT_BITS - mod.norm;

    /* at most k*min(amax, bmax) products, each < 2^(2*modbits), accumulate
       onto one output coefficient */
    {
        ulong lo_prod, hi_prod;
        umul_ppmm(hi_prod, lo_prod, (ulong) k, FLINT_MIN(amax, bmax));
        if (hi_prod != 0)
            return 0;
        len_bound = lo_prod;
    }

    success = fft_small_plan_init_nmod(P, get_default_mpn_ctx(),
                    (ulong) zl, hi, zn,
                    n_max(n_round_up(amax, BLK_SZ), n_round_up(bmax, BLK_SZ)),
                    len_bound, 2*modbits, mod, FLINT_MIN(amax, bmax));
    if (!success)
        return 0;

    nA = (ulong) ar * k;
    nB = (ulong) k * bc;
    nC = (ulong) ar * bc;

    ops = FLINT_ARRAY_ALLOC(nA + nB, fft_small_op_struct);
    tasks = FLINT_ARRAY_ALLOC(nA + nB, const nmod_poly_struct*);

    for (i = 0; i < ar; i++)
        for (l = 0; l < k; l++)
            tasks[i*k + l] = nmod_poly_mat_entry(A, i, l);
    for (l = 0; l < k; l++)
        for (j = 0; j < bc; j++)
            tasks[nA + l*bc + j] = nmod_poly_mat_entry(B, l, j);

    /* itrunc == 0 marks an untransformed (zero) operand for the
       accumulation pass */
    for (t = 0; t < nA + nB; t++)
    {
        fft_small_op_init(ops + t, P);
        ops[t].itrunc = 0;
    }

    ntasks = FLINT_MAX(nA + nB, nC);
    nworkers = flint_request_threads(&handles, ntasks);
    nthreads = nworkers + 1;

    /* stage 1: all input transforms, distributed over the threads; when
       tasks are fewer than threads, fft_small_fft_nmod itself picks up
       the spare pool threads */
    targs = FLINT_ARRAY_ALLOC(nthreads, _transform_worker_struct);
    for (t = 0; t < nthreads; t++)
    {
        _transform_worker_struct* W = targs + t;
        W->ops = ops;
        W->tasks = tasks;
        W->ntasks = nA + nB;
        W->mod = mod;
        W->P = P;
        W->start = (t+0)*(nA + nB)/nthreads;
        W->stop  = (t+1)*(nA + nB)/nthreads;
    }

    for (t = nthreads - 1; t > 0; t--)
        thread_pool_wake(global_thread_pool, handles[t - 1], 0,
                         _transform_worker_func, targs + t);
    _transform_worker_func(targs + 0);
    for (t = nthreads - 1; t > 0; t--)
        thread_pool_wait(global_thread_pool, handles[t - 1]);

    /* stage 2: accumulate, invert and export each output entry,
       distributed over the threads */
    eargs = FLINT_ARRAY_ALLOC(nthreads, _entry_worker_struct);
    for (t = 0; t < nthreads; t++)
    {
        _entry_worker_struct* W = eargs + t;
        W->C = C;
        W->Aops = ops;
        W->Bops = ops + nA;
        W->k = k;
        W->bc = bc;
        W->wlen = wlen;
        W->mod = mod;
        W->P = P;
        fft_small_op_init(W->Z, P);
        W->tmp = FLINT_ARRAY_ALLOC(P->zh - P->zl, ulong);
        W->start_ij = (t+0)*nC/nthreads;
        W->stop_ij  = (t+1)*nC/nthreads;
    }

    for (t = nthreads - 1; t > 0; t--)
        thread_pool_wake(global_thread_pool, handles[t - 1], 0,
                         _entry_worker_func, eargs + t);
    _entry_worker_func(eargs + 0);
    for (t = nthreads - 1; t > 0; t--)
        thread_pool_wait(global_thread_pool, handles[t - 1]);

    flint_give_back_threads(handles, nworkers);

    for (t = 0; t < nthreads; t++)
    {
        fft_small_op_clear(eargs[t].Z);
        flint_free(eargs[t].tmp);
    }
    for (t = 0; t < nA + nB; t++)
        fft_small_op_clear(ops + t);

    flint_free(eargs);
    flint_free(targs);
    flint_free(tasks);
    flint_free(ops);
    fft_small_plan_clear(P);

    return 1;
}

#else /* FLINT_HAVE_FFT_SMALL */

int nmod_poly_mat_mulmid_fft_small(nmod_poly_mat_t C,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                    slong zl, slong zh)
{
    return 0;
}

#endif
