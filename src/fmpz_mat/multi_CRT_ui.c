/*
    Copyright (C) 2010, 2011, 2018, 2026 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "mpn_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

typedef struct {
    slong n;
    slong Cstartrow;
    slong Cstoprow;
    fmpz * Centries;
    slong Cstride;
    nmod_mat_t * mod_C;
    nn_ptr primes;
    slong num_primes;
    const fmpz_comb_struct * comb;
    int sign;
} _crt_worker_arg;

static void _crt_worker(void * varg)
{
    _crt_worker_arg * arg = (_crt_worker_arg *) varg;
    slong i, j, l;
    slong n            = arg->n;
    slong Cstartrow    = arg->Cstartrow;
    slong Cstoprow     = arg->Cstoprow;
    fmpz * Centries    = arg->Centries;
    slong Cstride      = arg->Cstride;
    nmod_mat_t * mod_C = arg->mod_C;
    slong num_primes   = arg->num_primes;
    const fmpz_comb_struct * comb = arg->comb;
    int sign           = arg->sign;

    FLINT_ASSERT(sign == 0 || sign == 1);

    if (comb != NULL)
    {
        ulong * residues;
        fmpz_comb_temp_t comb_temp;

        residues = FLINT_ARRAY_ALLOC(num_primes, ulong);
        fmpz_comb_temp_init(comb_temp, comb);

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            for (l = 0; l < num_primes; l++)
                residues[l] = nmod_mat_entry(mod_C[l], i, j);
            fmpz_multi_CRT_ui(Centries + i * Cstride + j, residues, comb, comb_temp, sign);
        }

        flint_free(residues);
        fmpz_comb_temp_clear(comb_temp);
    }
    else if (num_primes == 1)
    {
        ulong r, t, p = mod_C[0]->mod.n;

        if (sign)
        {
            for (i = Cstartrow; i < Cstoprow; i++)
            for (j = 0; j < n; j++)
            {
                r = nmod_mat_entry(mod_C[0], i, j);
                t = p - r;
                if (t < r)
                    fmpz_neg_ui(Centries + i * Cstride + j, t);
                else
                    fmpz_set_ui(Centries + i * Cstride + j, r);
            }
        }
        else
        {
            for (i = Cstartrow; i < Cstoprow; i++)
            for (j = 0; j < n; j++)
            {
                r = nmod_mat_entry(mod_C[0], i, j);
                fmpz_set_ui(Centries + i * Cstride + j, r);
            }
        }
    }
    else if (num_primes == 2)
    {
        ulong r0, r1, i0, i1, m0, m1, M[2], t[2], u[2];
        m0 = mod_C[0]->mod.n;
        m1 = mod_C[1]->mod.n;
        umul_ppmm(M[1], M[0], m0, m1);
        if (FLINT_BIT_COUNT(M[1]) >= FLINT_BITS)
            goto generic;

        i0 = n_invmod(m1 % m0, m0);
        i1 = n_invmod(m0 % m1, m1);

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            r0 = nmod_mul(i0, nmod_mat_entry(mod_C[0], i, j), mod_C[0]->mod);
            r1 = nmod_mul(i1, nmod_mat_entry(mod_C[1], i, j), mod_C[1]->mod);

            FLINT_ASSERT(FLINT_BIT_COUNT(M[1]) < FLINT_BITS);
            umul_ppmm(t[1], t[0], r0, m1);
            umul_ppmm(u[1], u[0], r1, m0);
            add_ssaaaa(t[1], t[0], t[1], t[0], u[1], u[0]);

            if (t[1] > M[1] || (t[1] == M[1] && t[0] > M[0]))
                sub_ddmmss(t[1], t[0], t[1], t[0], M[1], M[0]);

            FLINT_ASSERT(t[1] < M[1] || (t[1] == M[1] && t[0] < M[0]));

            if (sign)
            {
                sub_ddmmss(u[1], u[0], M[1], M[0], t[1], t[0]);
                if (u[1] < t[1] || (u[1] == t[1] && u[0] < t[0]))
                    fmpz_neg_uiui(Centries + i * Cstride + j, u[1], u[0]);
                else
                    fmpz_set_uiui(Centries + i * Cstride + j, t[1], t[0]);
            }
            else
            {
                fmpz_set_uiui(Centries + i * Cstride + j, t[1], t[0]);
            }
        }
    }
    else
    {
        nn_ptr M, Ns, T, U, primes;
        slong Msize, Nsize;
        ulong cy, ri;

generic:
        primes = FLINT_ARRAY_ALLOC(num_primes, ulong);
        M = FLINT_ARRAY_ALLOC(num_primes + 1, ulong);
        M[0] = primes[0] = mod_C[0]->mod.n;
        Msize = 1;
        for (i = 1; i < num_primes; i++)
        {
            primes[i] = mod_C[i]->mod.n;
            M[Msize] = cy = mpn_mul_1(M, M, Msize, primes[i]);
            Msize += (cy != 0);
        }

        Nsize = Msize + 2;
        Ns = FLINT_ARRAY_ALLOC(Nsize * num_primes, ulong);
        T  = FLINT_ARRAY_ALLOC(Nsize, ulong);
        U  = FLINT_ARRAY_ALLOC(Nsize, ulong);

        for (i = 0; i < num_primes; i++)
        {
            Ns[i*Nsize + (Nsize - 1)] = 0;
            Ns[i*Nsize + (Nsize - 2)] = 0;
            mpn_divexact_1(Ns + i * Nsize, M, Msize, primes[i]);
            ri = mpn_mod_1(Ns + i * Nsize, Msize, primes[i]);
            ri = n_invmod(ri, primes[i]);
            Ns[i*Nsize + Msize] = mpn_mul_1(Ns + i*Nsize, Ns + i*Nsize, Msize, ri);
        }

        for (i = Cstartrow; i < Cstoprow; i++)
        for (j = 0; j < n; j++)
        {
            ri = nmod_mat_entry(mod_C[0], i, j);
            T[Nsize - 1] = mpn_mul_1(T, Ns, Nsize - 1, ri);
            for (l = 1; l < num_primes; l++)
            {
                ri = nmod_mat_entry(mod_C[l], i, j);
                T[Nsize - 1] += mpn_addmul_1(T, Ns + l*Nsize, Nsize - 1, ri);
            }
            mpn_tdiv_qr(U, T, 0, T, Nsize, M, Msize);
            if (sign && (mpn_sub_n(U, M, T, Msize), mpn_cmp(U, T, Msize) < 0))
            {
                fmpz_set_ui_array(Centries + i * Cstride + j, U, Msize);
                fmpz_neg(Centries + i * Cstride + j, Centries + i * Cstride + j);
            }
            else
            {
                fmpz_set_ui_array(Centries + i * Cstride + j, T, Msize);
            }
        }

        flint_free(primes);
        flint_free(M);
        flint_free(Ns);
        flint_free(T);
        flint_free(U);
    }
}

void fmpz_mat_multi_CRT_ui_precomp(
    fmpz_mat_t mat,
    nmod_mat_t * const residues, slong nres,
    const fmpz_comb_t comb,
    int sign)
{
    slong i, start, stop;
    slong m = fmpz_mat_nrows(mat);
    slong n = fmpz_mat_ncols(mat);
    slong num_workers, thread_limit;
    thread_pool_handle * handles;
    _crt_worker_arg mainarg, * args;

    mainarg.n          = n;
    mainarg.Centries   = mat->entries;
    mainarg.Cstride    = mat->stride;
    mainarg.mod_C      = residues;
    mainarg.num_primes = nres;
    mainarg.comb       = comb;
    mainarg.sign       = sign;

    thread_limit = (m * n * nres) / FMPZ_MAT_CRT_MIN_WORK_PER_THREAD;

    if (thread_limit < 2)
    {
single:
        mainarg.Cstartrow = 0;
        mainarg.Cstoprow  = m;
        _crt_worker(&mainarg);
        return;
    }

    num_workers = flint_request_threads(&handles, thread_limit);
    if (num_workers < 1)
    {
        flint_give_back_threads(handles, num_workers);
        goto single;
    }

    args = FLINT_ARRAY_ALLOC(num_workers, _crt_worker_arg);
    for (start = 0, i = 0; i < num_workers; start = stop, i++)
    {
        args[i]           = mainarg;
        stop              = (i + 1) * m / (num_workers + 1);
        args[i].Cstartrow = start;
        args[i].Cstoprow  = stop;
    }
    mainarg.Cstartrow = start;
    mainarg.Cstoprow  = m;

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _crt_worker, &args[i]);
    _crt_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);
}

void fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, int sign)
{
    if (nres > FMPZ_MAT_CRT_PRIMES_COMB_CUTOFF)
    {
        fmpz_comb_t comb;
        slong i;
        nn_ptr primes;
        primes = _nmod_vec_init(nres);
        for (i = 0; i < nres; i++)
            primes[i] = residues[i]->mod.n;
        fmpz_comb_init(comb, primes, nres);
        _nmod_vec_clear(primes);

        fmpz_mat_multi_CRT_ui_precomp(mat, residues, nres, comb, sign);

        fmpz_comb_clear(comb);
    }
    else
    {
        fmpz_mat_multi_CRT_ui_precomp(mat, residues, nres, NULL, sign);
    }
}

