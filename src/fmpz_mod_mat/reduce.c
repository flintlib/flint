/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"


typedef struct {
    slong startrow;
    slong stoprow;
    fmpz_mod_mat_struct * M;
} _worker_arg;

static void _red_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong startrow = arg->startrow;
    slong stoprow = arg->stoprow;
    fmpz_mod_mat_struct * M = arg->M;
    slong c = fmpz_mod_mat_ncols(M);
    slong i;

    for (i = startrow; i < stoprow; i++)
        _fmpz_vec_scalar_mod_fmpz(M->mat->rows[i], M->mat->rows[i], c, M->mod);
}

void _fmpz_mod_mat_reduce(fmpz_mod_mat_t M)
{
    slong i, r = fmpz_mod_mat_nrows(M);
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg mainarg;
    _worker_arg * args;
    slong limit;

    /* limit on threads */
    limit = fmpz_size(M->mod) + r + fmpz_mod_mat_ncols(M);
    limit = limit < 64 ? 0 : (limit - 64)/64;
    limit = FLINT_MIN(limit, r);

    mainarg.startrow = 0;
    mainarg.stoprow = r;
    mainarg.M = M;

    if (limit < 2)
    {
use_one_thread:
        _red_worker(&mainarg);
        return;
    }

    num_workers = flint_request_threads(&handles, limit);
    if (num_workers < 1)
    {
        flint_give_back_threads(handles, num_workers);
        goto use_one_thread;
    }

    args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);

    for (i = 0; i < num_workers; i++)
    {
        args[i].startrow = (i + 0)*r/(num_workers + 1);
        args[i].stoprow = (i + 1)*r/(num_workers + 1);
        args[i].M = M;
    }

    i = num_workers;
    mainarg.startrow = (i + 0)*r/(num_workers + 1);
    mainarg.stoprow = (i + 1)*r/(num_workers + 1);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _red_worker, &args[i]);
    _red_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);

    return;
}

