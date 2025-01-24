/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "acb.h"
#include "acb_mat.h"

typedef struct
{
    acb_ptr C;
    acb_srcptr A;
    acb_srcptr B;
    slong Cstride;
    slong Astride;
    slong Bstride;
    slong ar0;
    slong ar1;
    slong bc0;
    slong bc1;
    slong br;
    slong prec;
}
_worker_arg;

void
_acb_mat_mul_thread(void * arg_ptr)
{
    _worker_arg arg = *((_worker_arg *) arg_ptr);
    slong i, j, br;

    br = arg.br;

    for (i = arg.ar0; i < arg.ar1; i++)
    {
        for (j = arg.bc0; j < arg.bc1; j++)
        {
            acb_dot(arg.C + i * arg.Cstride + j, NULL, 0,
                arg.A + i * arg.Astride, 1,
                arg.B + j, arg.Bstride, br, arg.prec);
        }
    }

    flint_cleanup();
    return;
}

void
acb_mat_mul_threaded(acb_mat_t C, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong ar, ac, br, bc, i, num_threads, num_workers;
    thread_pool_handle * handles;
    _worker_arg * args;

    ar = acb_mat_nrows(A);
    ac = acb_mat_ncols(A);
    br = acb_mat_nrows(B);
    bc = acb_mat_ncols(B);

    if (ac != br || ar != acb_mat_nrows(C) || bc != acb_mat_ncols(C))
        flint_throw(FLINT_DOMERR, "incompatible dimensions in %s\n", __func__);

    if (br == 0)
    {
        acb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        acb_mat_t T;
        acb_mat_init(T, ar, bc);
        acb_mat_mul_threaded(T, A, B, prec);
        acb_mat_swap_entrywise(T, C);
        acb_mat_clear(T);
        return;
    }

    num_workers = flint_request_threads(&handles, FLINT_MAX(ar, bc));
    num_threads = num_workers + 1;

    args = FLINT_ARRAY_ALLOC(num_threads, _worker_arg);

    for (i = 0; i < num_threads; i++)
    {
        args[i].C = C->entries;
        args[i].A = A->entries;
        args[i].B = B->entries;
        args[i].Cstride = C->stride;
        args[i].Astride = A->stride;
        args[i].Bstride = B->stride;

        if (ar >= bc)
        {
            args[i].ar0 = (ar * i) / num_threads;
            args[i].ar1 = (ar * (i + 1)) / num_threads;
            args[i].bc0 = 0;
            args[i].bc1 = bc;
        }
        else
        {
            args[i].ar0 = 0;
            args[i].ar1 = ar;
            args[i].bc0 = (bc * i) / num_threads;
            args[i].bc1 = (bc * (i + 1)) / num_threads;
        }

        args[i].br = br;
        args[i].prec = prec;

        if (i < num_workers)
            thread_pool_wake(global_thread_pool, handles[i], 0, _acb_mat_mul_thread, &args[i]);
        else
            _acb_mat_mul_thread(&args[i]);
    }

    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);
}
