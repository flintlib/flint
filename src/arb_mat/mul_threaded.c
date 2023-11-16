/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "arb_mat.h"

typedef struct
{
    arb_ptr * C;
    const arb_ptr * A;
    const arb_ptr * B;
    slong ar0;
    slong ar1;
    slong bc0;
    slong bc1;
    slong br;
    slong prec;
}
_worker_arg;

void
_arb_mat_mul_thread(void * arg_ptr)
{
    _worker_arg arg = *((_worker_arg *) arg_ptr);
    slong i, j, br, bc;
    arb_ptr tmp;
    TMP_INIT;

    br = arg.br;
    bc = arg.bc1 - arg.bc0;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(arb_struct) * br * bc);

    for (i = 0; i < br; i++)
        for (j = 0; j < bc; j++)
            tmp[j * br + i] = arg.B[i][arg.bc0 + j];

    for (i = arg.ar0; i < arg.ar1; i++)
    {
        for (j = arg.bc0; j < arg.bc1; j++)
        {
            arb_dot(arg.C[i] + j, NULL, 0,
                arg.A[i], 1, tmp + (j - arg.bc0) * br, 1, br, arg.prec);
        }
    }

    TMP_END;
    flint_cleanup();
    return;
}

void
arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    slong ar, ac, br, bc, i, num_threads, num_workers;
    thread_pool_handle * handles;
    _worker_arg * args;

    ar = arb_mat_nrows(A);
    ac = arb_mat_ncols(A);
    br = arb_mat_nrows(B);
    bc = arb_mat_ncols(B);

    if (ac != br || ar != arb_mat_nrows(C) || bc != arb_mat_ncols(C))
        flint_throw(FLINT_DOMERR, "incompatible dimensions in %s\n", __func__);

    if (br == 0)
    {
        arb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        arb_mat_t T;
        arb_mat_init(T, ar, bc);
        arb_mat_mul_threaded(T, A, B, prec);
        arb_mat_swap_entrywise(T, C);
        arb_mat_clear(T);
        return;
    }

    num_workers = flint_request_threads(&handles, FLINT_MAX(ar, bc));
    num_threads = num_workers + 1;

    args = FLINT_ARRAY_ALLOC(num_threads, _worker_arg);

    for (i = 0; i < num_threads; i++)
    {
        args[i].C = C->rows;
        args[i].A = A->rows;
        args[i].B = B->rows;

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
            thread_pool_wake(global_thread_pool, handles[i], 0, _arb_mat_mul_thread, &args[i]);
        else
            _arb_mat_mul_thread(&args[i]);
    }

    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);
}
