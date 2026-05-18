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
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "mpn_extras.h"

typedef struct {
    slong ncols_A;
    slong ncols_B;
    slong Astartrow;
    slong Astoprow;
    slong Bstartrow;
    slong Bstoprow;
    fmpz * Aentries;
    slong Astride;
    fmpz * Bentries;
    slong Bstride;
    nmod_mat_t * mod_A;
    nmod_mat_t * mod_B;   /* NULL for single-matrix case */
    const fmpz_comb_struct * comb;
    slong num_primes;
} _mod_worker_arg;

static void _mod_worker(void * varg)
{
    _mod_worker_arg * arg = (_mod_worker_arg *) varg;
    slong i, j, l;
    slong ncols_A    = arg->ncols_A;
    slong ncols_B    = arg->ncols_B;
    slong Astartrow  = arg->Astartrow;
    slong Astoprow   = arg->Astoprow;
    slong Bstartrow  = arg->Bstartrow;
    slong Bstoprow   = arg->Bstoprow;
    fmpz * Aentries  = arg->Aentries;
    slong Astride    = arg->Astride;
    fmpz * Bentries  = arg->Bentries;
    slong Bstride    = arg->Bstride;
    nmod_mat_t * mod_A = arg->mod_A;
    nmod_mat_t * mod_B = arg->mod_B;
    slong num_primes   = arg->num_primes;
    const fmpz_comb_struct * comb = arg->comb;

    if (comb != NULL)
    {
        ulong * residues;
        fmpz_comb_temp_t comb_temp;

        residues = FLINT_ARRAY_ALLOC(num_primes, ulong);
        fmpz_comb_temp_init(comb_temp, comb);

        for (i = Astartrow; i < Astoprow; i++)
            for (j = 0; j < ncols_A; j++)
            {
                fmpz_multi_mod_ui(residues, Aentries + i * Astride + j, comb, comb_temp);
                for (l = 0; l < num_primes; l++)
                    nmod_mat_entry(mod_A[l], i, j) = residues[l];
            }

        if (mod_B != NULL)
            for (i = Bstartrow; i < Bstoprow; i++)
                for (j = 0; j < ncols_B; j++)
                {
                    fmpz_multi_mod_ui(residues, Bentries + i * Bstride + j, comb, comb_temp);
                    for (l = 0; l < num_primes; l++)
                        nmod_mat_entry(mod_B[l], i, j) = residues[l];
                }

        flint_free(residues);
        fmpz_comb_temp_clear(comb_temp);
    }
    else
    {
        for (i = Astartrow; i < Astoprow; i++)
            for (j = 0; j < ncols_A; j++)
                for (l = 0; l < num_primes; l++)
                    nmod_mat_entry(mod_A[l], i, j) = fmpz_get_nmod(
                        Aentries + i * Astride + j, mod_A[l]->mod);

        if (mod_B != NULL)
            for (i = Bstartrow; i < Bstoprow; i++)
                for (j = 0; j < ncols_B; j++)
                    for (l = 0; l < num_primes; l++)
                        nmod_mat_entry(mod_B[l], i, j) = fmpz_get_nmod(
                            Bentries + i * Bstride + j, mod_B[l]->mod);
    }
}

static void _fmpz_mat_multi_mod_ui_internal(
    nmod_mat_t * residues,
    nmod_mat_t * residues_B,   /* NULL for single-matrix case */
    const fmpz_mat_t A,
    const fmpz_mat_t B,        /* NULL for single-matrix case */
    const fmpz_comb_t comb,
    slong num_primes)
{
    slong i, start, stop;
    slong m  = A->r;
    slong k  = A->c;
    slong kB = (B != NULL) ? B->c : 0;
    slong num_workers, thread_limit;
    thread_pool_handle * handles;
    _mod_worker_arg mainarg, * args;

    mainarg.ncols_A   = k;
    mainarg.ncols_B   = kB;
    mainarg.Aentries  = A->entries;
    mainarg.Astride   = A->stride;
    mainarg.Bentries  = (B != NULL) ? B->entries : NULL;
    mainarg.Bstride   = (B != NULL) ? B->stride  : 0;
    mainarg.mod_A     = residues;
    mainarg.mod_B     = residues_B;
    mainarg.comb      = comb;
    mainarg.num_primes = num_primes;

    thread_limit = (FLINT_MAX(m * k, k * kB) * num_primes) / FMPZ_MAT_CRT_MIN_WORK_PER_THREAD;

    if (thread_limit < 2)
    {
single:
        mainarg.Astartrow = 0;
        mainarg.Astoprow  = m;
        mainarg.Bstartrow = 0;
        mainarg.Bstoprow  = (B != NULL) ? B->r : 0;
        _mod_worker(&mainarg);
        return;
    }

    num_workers = flint_request_threads(&handles, thread_limit);
    if (num_workers < 1)
    {
        flint_give_back_threads(handles, num_workers);
        goto single;
    }

    args = FLINT_ARRAY_ALLOC(num_workers, _mod_worker_arg);
    for (start = 0, i = 0; i < num_workers; start = stop, i++)
    {
        args[i] = mainarg;
        stop = _thread_pool_find_work_2(m, k, (B != NULL) ? B->r : 0, kB,
                                        i + 1, num_workers + 1);
        _thread_pool_distribute_work_2(start, stop,
                                       &args[i].Astartrow, &args[i].Astoprow, m,
                                       &args[i].Bstartrow, &args[i].Bstoprow,
                                       (B != NULL) ? B->r : 0);
    }
    _thread_pool_distribute_work_2(start, m + ((B != NULL) ? B->r : 0),
                                   &mainarg.Astartrow, &mainarg.Astoprow, m,
                                   &mainarg.Bstartrow, &mainarg.Bstoprow,
                                   (B != NULL) ? B->r : 0);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _mod_worker, &args[i]);
    _mod_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);
}

void fmpz_mat_multi_mod_ui_precomp(
    nmod_mat_t * residues, slong nres,
    const fmpz_mat_t A,
    const fmpz_comb_t comb)
{
    _fmpz_mat_multi_mod_ui_internal(residues, NULL, A, NULL, comb, nres);
}

void fmpz_mat_multi_mod_2_ui_precomp(
    nmod_mat_t * residues_A, nmod_mat_t * residues_B, slong nres,
    const fmpz_mat_t A, const fmpz_mat_t B,
    const fmpz_comb_t comb)
{
    _fmpz_mat_multi_mod_ui_internal(residues_A, residues_B, A, B, comb, nres);
}

void fmpz_mat_multi_mod_ui(nmod_mat_t * residues, slong nres, const fmpz_mat_t A)
{
    if (nres > FMPZ_MAT_MOD_PRIMES_COMB_CUTOFF)
    {
        fmpz_comb_t comb;
        slong i;
        nn_ptr primes;
        primes = _nmod_vec_init(nres);
        for (i = 0; i < nres; i++)
            primes[i] = residues[i]->mod.n;
        fmpz_comb_init(comb, primes, nres);
        _nmod_vec_clear(primes);

        fmpz_mat_multi_mod_ui_precomp(residues, nres, A, comb);

        fmpz_comb_clear(comb);
    }
    else
    {
        fmpz_mat_multi_mod_ui_precomp(residues, nres, A, NULL);
    }
}

