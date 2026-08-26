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
    fmpz_comb_t comb_local;
    int own_comb = 0;

    if (comb == NULL)
    {
        nn_ptr primes = FLINT_ARRAY_ALLOC(nres, ulong);
        for (i = 0; i < nres; i++)
            primes[i] = residues[i]->mod.n;
        fmpz_comb_init2(comb_local, primes, nres, FMPZ_COMB_CRT);
        flint_free(primes);
        comb = comb_local;
        own_comb = 1;
    }

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
        goto cleanup;
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

cleanup:
    if (own_comb)
        fmpz_comb_clear(comb_local);
}

void fmpz_mat_multi_CRT_ui(fmpz_mat_t mat, nmod_mat_t * const residues, slong nres, int sign)
{
    fmpz_mat_multi_CRT_ui_precomp(mat, residues, nres, NULL, sign);
}
