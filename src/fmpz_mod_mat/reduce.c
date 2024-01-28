/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"

/* todo: rewrite using parallel_do */

typedef struct {
    slong startrow;
    slong stoprow;
    fmpz_mod_mat_struct * M;
    const fmpz_mod_ctx_struct * ctx;
} _worker_arg;

static void _red_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong startrow = arg->startrow;
    slong stoprow = arg->stoprow;
    fmpz_mod_mat_struct * M = arg->M;
    const fmpz_mod_ctx_struct * ctx = arg->ctx;
    slong c = fmpz_mod_mat_ncols(M, ctx);
    slong i;

    for (i = startrow; i < stoprow; i++)
        _fmpz_mod_vec_set_fmpz_vec(M->rows[i], M->rows[i], c, ctx);
}

void _fmpz_mod_mat_reduce(fmpz_mod_mat_t M, const fmpz_mod_ctx_t ctx)
{
    slong i, r = fmpz_mod_mat_nrows(M, ctx);
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg mainarg;
    _worker_arg * args;
    slong limit;

    /* limit on threads */
    limit = fmpz_size(ctx->n) + r + fmpz_mod_mat_ncols(M, ctx);
    limit = limit < 64 ? 0 : (limit - 64)/64;
    limit = FLINT_MIN(limit, r);

    mainarg.startrow = 0;
    mainarg.stoprow = r;
    mainarg.M = M;
    mainarg.ctx = ctx;

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
        args[i].ctx = ctx;
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
