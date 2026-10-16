/*
    Copyright (C) 2024 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Multimodular matrix multiplication largely copied from fmpz_mat.
   Should think about how to reuse code here and with other types. */

#include "thread_pool.h"
#include "thread_support.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "mpn_extras.h"
#include "mpn_mod.h"
#include "gr_mat.h"

typedef struct {
    slong m;
    slong k;
    slong n;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    slong Cstartrow;
    slong Cstoprow;
    nn_ptr Aentries;
    slong Astride;
    nn_ptr Bentries;
    slong Bstride;
    nn_ptr Centries;
    slong Cstride;
    nmod_mat_t * mod_A;
    nmod_mat_t * mod_B;
    nmod_mat_t * mod_C;
    slong num_primes;
    nn_ptr primes;
    const flint_mpn_crt_struct * crt;
    gr_ctx_struct * ctx;
} _worker_arg;

static void _mod_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong i, j, l;
    slong k = arg->k;
    slong n = arg->n;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong Bstartrow = arg->Bstartrow;
    slong Bstoprow = arg->Bstoprow;
    nn_ptr Aentries = arg->Aentries;
    slong Astride = arg->Astride;
    nn_ptr Bentries = arg->Bentries;
    slong Bstride = arg->Bstride;
    nmod_mat_t * mod_A = arg->mod_A;
    nmod_mat_t * mod_B = arg->mod_B;
    slong num_primes = arg->num_primes;
    const flint_mpn_crt_struct * crt = arg->crt;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(arg->ctx);
    nn_ptr res, tmp;
    slong maxcols = FLINT_MAX(k, (mod_B != NULL) ? n : 0);

    res = FLINT_ARRAY_ALLOC(num_primes * maxcols, ulong);
    tmp = FLINT_ARRAY_ALLOC(crt->tmp_limbs, ulong);

    /* the entries of a row are contiguous, so use the vector function */
    for (i = Astartrow; i < Astoprow; i++)
    {
        flint_mpn_multi_mod_vec(res, k, Aentries + i * Astride * nlimbs, nlimbs, k, crt, tmp);
        for (l = 0; l < num_primes; l++)
            for (j = 0; j < k; j++)
                nmod_mat_entry(mod_A[l], i, j) = res[l * k + j];
    }

    if (mod_B != NULL)
    {
        for (i = Bstartrow; i < Bstoprow; i++)
        {
            flint_mpn_multi_mod_vec(res, n, Bentries + i * Bstride * nlimbs, nlimbs, n, crt, tmp);
            for (l = 0; l < num_primes; l++)
                for (j = 0; j < n; j++)
                    nmod_mat_entry(mod_B[l], i, j) = res[l * n + j];
        }
    }

    flint_free(res);
    flint_free(tmp);
}

static void _crt_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong i, j, l;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    nn_ptr Centries = arg->Centries;
    slong Cstride = arg->Cstride;
    nmod_mat_t * mod_C = arg->mod_C;
    slong num_primes = arg->num_primes;
    const flint_mpn_crt_struct * crt = arg->crt;
    gr_ctx_struct * ctx = arg->ctx;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    nn_ptr res, out, tmp;
    slong pn = crt->prod_len;

    res = FLINT_ARRAY_ALLOC(num_primes * n, ulong);
    out = FLINT_ARRAY_ALLOC(n * pn, ulong);
    tmp = FLINT_ARRAY_ALLOC(crt->tmp_limbs, ulong);

    for (i = Cstartrow; i < Cstoprow; i++)
    {
        for (l = 0; l < num_primes; l++)
            for (j = 0; j < n; j++)
                res[l * n + j] = nmod_mat_entry(mod_C[l], i, j);

        flint_mpn_multi_crt_vec(out, pn, NULL, res, n, n, crt, 0, tmp);

        for (j = 0; j < n; j++)
            mpn_mod_set_mpn(Centries + (i * Cstride + j) * nlimbs, out + j * pn, pn, ctx);
    }

    flint_free(res);
    flint_free(out);
    flint_free(tmp);
}

int mpn_mod_mat_mul_multi_mod(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong i, start, stop;
    slong m, k, n;
    flint_bitcnt_t primes_bits;
    _worker_arg mainarg;
    _worker_arg * args;
    slong num_workers;
    thread_pool_handle * handles;
    slong limit;
    ulong first_prime; /* not prime */
    int squaring = (A == B);
    flint_bitcnt_t Abits, Bbits, Cbits, bits, mod_bits;
    slong nlimbs;

    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    mod_bits = FLINT_BITS * (nlimbs - 1) + (FLINT_BITS - MPN_MOD_CTX_NORM(ctx));

    /* todo: consider optimizing for small entries */
    Abits = mod_bits;
    Bbits = mod_bits;
    Cbits = Abits + Bbits + FLINT_BIT_COUNT(A->c);
    bits = Cbits;

    mainarg.m = m = A->r;
    mainarg.k = k = A->c;
    mainarg.n = n = B->c;

    if (m < 1 || n < 1 || k < 1)
        return gr_mat_zero(C, ctx);

    mainarg.ctx = ctx;

    mainarg.Aentries = (nn_ptr) A->entries;
    mainarg.Astride = A->stride;
    mainarg.Bentries = (nn_ptr) B->entries;
    mainarg.Bstride = B->stride;
    mainarg.Centries = (nn_ptr) C->entries;
    mainarg.Cstride = C->stride;

    /* TUNING */
    primes_bits = NMOD_MAT_OPTIMAL_MODULUS_BITS;

    if (bits < primes_bits || bits <= FLINT_BITS - 1)
    {
        mainarg.num_primes = 1;
        first_prime = UWORD(1) << bits;
    }
    else
    {
        /* Round up in the division. The first modulus is a power of two
           of the same size as the primes, so that the fixed-length CRT
           code in flint_mpn_multi_crt applies for small products. */
        mainarg.num_primes = 1 + (bits - primes_bits + primes_bits - 1)/primes_bits;
        first_prime = UWORD(1) << primes_bits;
    }

    /* Initialize */
    mainarg.primes = FLINT_ARRAY_ALLOC(mainarg.num_primes, ulong);
    mainarg.primes[0] = first_prime;
    if (mainarg.num_primes > 1)
    {
        mainarg.primes[1] = n_nextprime(UWORD(1) << primes_bits, 0);
        for (i = 2; i < mainarg.num_primes; i++)
            mainarg.primes[i] = n_nextprime(mainarg.primes[i-1], 0);
    }

    {
        flint_mpn_crt_struct * crt = flint_malloc(sizeof(flint_mpn_crt_struct));
        flint_mpn_crt_init(crt, mainarg.primes, mainarg.num_primes);
        mainarg.crt = crt;
    }

    mainarg.mod_A = FLINT_ARRAY_ALLOC(mainarg.num_primes, nmod_mat_t);

    if (squaring)
        mainarg.mod_B = NULL;
    else
        mainarg.mod_B = FLINT_ARRAY_ALLOC(mainarg.num_primes, nmod_mat_t);

    mainarg.mod_C = FLINT_ARRAY_ALLOC(mainarg.num_primes, nmod_mat_t);
    for (i = 0; i < mainarg.num_primes; i++)
    {
        nmod_mat_init(mainarg.mod_A[i], A->r, A->c, mainarg.primes[i]);
        if (!squaring)
            nmod_mat_init(mainarg.mod_B[i], B->r, B->c, mainarg.primes[i]);
        nmod_mat_init(mainarg.mod_C[i], C->r, C->c, mainarg.primes[i]);
    }

    /* limit on the number of threads */
    limit = ((m + k + n)/128)*(1 + bits/1024);
    limit = FLINT_MIN(limit, (m + k)/4);

    /* mod */
    if (limit < 2)
    {
mod_single:
        mainarg.Astartrow = 0;
        mainarg.Astoprow = m;
        mainarg.Bstartrow = 0;
        mainarg.Bstoprow = k;
        _mod_worker(&mainarg);
    }
    else
    {
        num_workers = flint_request_threads(&handles, limit);
        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);
            goto mod_single;
        }

        args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);
        for (start = 0, i = 0; i < num_workers; start = stop, i++)
        {
            args[i] = mainarg;
            stop = _thread_pool_find_work_2(m, k, k, n, i + 1, num_workers + 1);
            _thread_pool_distribute_work_2(start, stop,
                                     &args[i].Astartrow, &args[i].Astoprow, m,
                                     &args[i].Bstartrow, &args[i].Bstoprow, k);
        }

        _thread_pool_distribute_work_2(start, m + k,
                                     &mainarg.Astartrow, &mainarg.Astoprow, m,
                                     &mainarg.Bstartrow, &mainarg.Bstoprow, k);

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _mod_worker, &args[i]);
        _mod_worker(&mainarg);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
        flint_free(args);
    }

    /* mul */
    for (i = 0; i < mainarg.num_primes; i++)
        nmod_mat_mul(mainarg.mod_C[i], mainarg.mod_A[i], squaring ? mainarg.mod_A[i] : mainarg.mod_B[i]);

    /* limit on the number of threads */
    limit = ((m + n)/64)*(1 + bits/1024);
    limit = FLINT_MIN(limit, m/2);

    /* crt */
    if (limit < 2)
    {
crt_single:
        mainarg.Cstartrow = 0;
        mainarg.Cstoprow = m;
        _crt_worker(&mainarg);
    }
    else
    {
        num_workers = flint_request_threads(&handles, limit);
        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);
            goto crt_single;
        }

        args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);
        for (start = 0, i = 0; i < num_workers; start = stop, i++)
        {
            args[i] = mainarg;
            stop = (i + 1)*m/(num_workers + 1);
            args[i].Cstartrow = start;
            args[i].Cstoprow = stop;
        }

        mainarg.Cstartrow = start;
        mainarg.Cstoprow = m;

        for (i = 0; i < num_workers; i++)
            thread_pool_wake(global_thread_pool, handles[i], 0, _crt_worker, &args[i]);
        _crt_worker(&mainarg);
        for (i = 0; i < num_workers; i++)
            thread_pool_wait(global_thread_pool, handles[i]);

        flint_give_back_threads(handles, num_workers);
        flint_free(args);
    }

    for (i = 0; i < mainarg.num_primes; i++)
    {
        nmod_mat_clear(mainarg.mod_A[i]);
        if (!squaring)
            nmod_mat_clear(mainarg.mod_B[i]);
        nmod_mat_clear(mainarg.mod_C[i]);
    }

    flint_free(mainarg.mod_A);
    if (!squaring)
        flint_free(mainarg.mod_B);
    flint_free(mainarg.mod_C);
    flint_free(mainarg.primes);
    flint_mpn_crt_clear((flint_mpn_crt_struct *) mainarg.crt);
    flint_free((flint_mpn_crt_struct *) mainarg.crt);

    return GR_SUCCESS;
}
