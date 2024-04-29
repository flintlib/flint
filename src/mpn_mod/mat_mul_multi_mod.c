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
#include "mpn_mod.h"

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
    mp_ptr * Arows;
    mp_ptr * Brows;
    mp_ptr * Crows;
    nmod_mat_t * mod_A;
    nmod_mat_t * mod_B;
    nmod_mat_t * mod_C;
    slong num_primes;
    mp_ptr primes;
    gr_ctx_struct * ctx;
} _worker_arg;

FLINT_FORCE_INLINE mp_limb_t
nmod_set_mpn_2(mp_srcptr ad, nmod_t mod)
{
    mp_limb_t r = 0;
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

FLINT_FORCE_INLINE mp_limb_t
nmod_set_mpn_3(mp_srcptr ad, nmod_t mod)
{
    mp_limb_t r = 0;
    NMOD_RED2(r, r, ad[2], mod);
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

FLINT_FORCE_INLINE mp_limb_t
nmod_set_mpn_4(mp_srcptr ad, nmod_t mod)
{
    mp_limb_t r = 0;
    NMOD_RED2(r, r, ad[3], mod);
    NMOD_RED2(r, r, ad[2], mod);
    NMOD_RED2(r, r, ad[1], mod);
    NMOD_RED2(r, r, ad[0], mod);
    return r;
}

/* todo: precomputed inverse */
FLINT_FORCE_INLINE mp_limb_t
nmod_set_mpn(mp_srcptr ad, mp_size_t an, nmod_t mod)
{
    return mpn_mod_1(ad, an, mod.n);
}

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
    mp_ptr * Arows = arg->Arows;
    mp_ptr * Brows = arg->Brows;
    nmod_mat_t * mod_A = arg->mod_A;
    nmod_mat_t * mod_B = arg->mod_B;
    slong num_primes = arg->num_primes;

    slong nlimbs = MPN_MOD_CTX_NLIMBS(arg->ctx);

    mp_limb_t first_prime = UWORD(1) << (FLINT_BITS - 1);

    if (nlimbs == 2 && arg->primes[0] == first_prime)
    {
        for (i = Astartrow; i < Astoprow; i++)
        {
            for (j = 0; j < k; j++)
            {
                nmod_mat_entry(mod_A[0], i, j) = (Arows[i] + j * nlimbs)[0] & (first_prime - 1);
                for (l = 1; l < num_primes; l++)
                    nmod_mat_entry(mod_A[l], i, j) = nmod_set_mpn_2(Arows[i] + j * nlimbs, mod_A[l]->mod);
            }
        }

        if (mod_B != NULL)
        {
            for (i = Bstartrow; i < Bstoprow; i++)
                for (j = 0; j < n; j++)
                {
                    nmod_mat_entry(mod_B[0], i, j) = (Brows[i] + j * nlimbs)[0] & (first_prime - 1);
                    for (l = 1; l < num_primes; l++)
                        nmod_mat_entry(mod_B[l], i, j) = nmod_set_mpn_2(Brows[i] + j * nlimbs, mod_A[l]->mod);
                }
        }
    }
    else
    {
        for (i = Astartrow; i < Astoprow; i++)
            for (j = 0; j < k; j++)
                for (l = 0; l < num_primes; l++)
                    nmod_mat_entry(mod_A[l], i, j) = nmod_set_mpn(Arows[i] + j * nlimbs, nlimbs, mod_A[l]->mod);

        if (mod_B != NULL)
        {
            for (i = Bstartrow; i < Bstoprow; i++)
                for (j = 0; j < n; j++)
                    for (l = 0; l < num_primes; l++)
                        nmod_mat_entry(mod_B[l], i, j) = nmod_set_mpn(Brows[i] + j * nlimbs, nlimbs, mod_A[l]->mod);
        }
    }
}

static void _crt_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong i, j, l;
    slong n = arg->n;
    slong Cstartrow = arg->Cstartrow;
    slong Cstoprow = arg->Cstoprow;
    mp_ptr * Crows = arg->Crows;
    nmod_mat_t * mod_C = arg->mod_C;
    mp_limb_t * primes = arg->primes;
    slong num_primes = arg->num_primes;
    gr_ctx_struct * ctx = arg->ctx;
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    /* todo: fmpz_mat has dedicated code for num_primes <= 2
             which we currently don't because we do not optimize
             for small entries */

    {
        mp_ptr M, Ns, T, U;
        mp_size_t Msize, Nsize;
        mp_limb_t cy, ri;

        M = FLINT_ARRAY_ALLOC(num_primes + 1, mp_limb_t);

        M[0] = primes[0];
        Msize = 1;
        for (i = 1; i < num_primes; i++)
        {
            FLINT_ASSERT(Msize > 0);
            M[Msize] = cy = mpn_mul_1(M, M, Msize, primes[i]);
            Msize += (cy != 0);
        }

        /* We add terms with Msize + 1 limbs, with one extra limb for the
           carry accumulation. todo: reduce Nsize by 1 when the carries
           do not require an extra limb. */
        Nsize = Msize + 2;

        Ns = FLINT_ARRAY_ALLOC(Nsize*num_primes, mp_limb_t);
        T = FLINT_ARRAY_ALLOC(Nsize, mp_limb_t);
        U = FLINT_ARRAY_ALLOC(Nsize, mp_limb_t);

        for (i = 0; i < num_primes; i++)
        {
            Ns[i*Nsize + (Nsize - 1)] = 0;
            Ns[i*Nsize + (Nsize - 2)] = 0;
            mpn_divrem_1(Ns + i * Nsize, 0, M, Msize, primes[i]);
            ri = mpn_mod_1(Ns + i * Nsize, Msize, primes[i]);
            ri = n_invmod(ri, primes[i]);
            FLINT_ASSERT(Msize > 0);
            Ns[i*Nsize + Msize] = mpn_mul_1(Ns + i*Nsize, Ns + i*Nsize, Msize, ri);
        }

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            ri = nmod_mat_entry(mod_C[0], i, j);
            FLINT_ASSERT(Nsize > 1);
            T[Nsize - 1] = mpn_mul_1(T, Ns, Nsize - 1, ri);

            /* todo: more fixed-length code */
            for (l = 1; l < num_primes; l++)
            {
                ri = nmod_mat_entry(mod_C[l], i, j);
                T[Nsize - 1] += mpn_addmul_1(T, Ns + l*Nsize, Nsize - 1, ri);
            }

            /* todo: division with precomputed inverse? */
            /* todo: see what fft_small does here */
            mpn_tdiv_qr(U, T, 0, T, Nsize, M, Msize);

            mpn_mod_set_mpn(Crows[i] + j * nlimbs, T, Msize, ctx);
        }

        flint_free(M);
        flint_free(Ns);
        flint_free(T);
        flint_free(U);
    }
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

    mainarg.Arows = (mp_ptr *) A->rows;
    mainarg.Brows = (mp_ptr *) B->rows;
    mainarg.Crows = (mp_ptr *) C->rows;

    /* TUNING */
    primes_bits = NMOD_MAT_OPTIMAL_MODULUS_BITS;

    if (bits < primes_bits || bits <= FLINT_BITS - 1)
    {
        mainarg.num_primes = 1;
        first_prime = UWORD(1) << bits;
    }
    else
    {
        /* Round up in the division */
        mainarg.num_primes = 1 + (bits - (FLINT_BITS - 1) + primes_bits - 1)/primes_bits;
        first_prime = UWORD(1) << (FLINT_BITS - 1);
    }

    /* Initialize */
    mainarg.primes = FLINT_ARRAY_ALLOC(mainarg.num_primes, mp_limb_t);
    mainarg.primes[0] = first_prime;
    if (mainarg.num_primes > 1)
    {
        mainarg.primes[1] = n_nextprime(UWORD(1) << primes_bits, 0);
        for (i = 2; i < mainarg.num_primes; i++)
            mainarg.primes[i] = n_nextprime(mainarg.primes[i-1], 0);
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

    return GR_SUCCESS;
}
