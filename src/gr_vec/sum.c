/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "gr_vec.h"


typedef struct
{
    gr_method_vec_reduce_op f;
    gr_srcptr vec;
    gr_ptr res;
    slong a;
    slong b;
    slong step;
    gr_ctx_struct * ctx;
    int status;
}
work_chunk_t;

static void
worker(void * _work)
{
    work_chunk_t work = *((work_chunk_t *) _work);
    int status = GR_SUCCESS;

    status |= work.f(work.res, GR_ENTRY(work.vec, work.a, work.ctx->sizeof_elem), work.b - work.a, work.ctx);
    work.status = status;
}

int _gr_vec_parallel_reduce(gr_ptr res, gr_method_vec_reduce_op basecase, gr_srcptr vec, slong n, gr_ctx_t ctx, int thread_limit, int flags)
{
    int status = GR_SUCCESS;

    if (thread_limit <= 0)
        thread_limit = flint_get_num_threads();

    thread_limit = FLINT_MIN(thread_limit, n / 2);

    if (thread_limit <= 1)
    {
        return basecase(res, vec, n, ctx);
    }
    else
    {
        slong i, num_threads, num_workers;
        thread_pool_handle * handles;

        num_workers = flint_request_threads(&handles, thread_limit);
        num_threads = num_workers + 1;

        if (flags & FLINT_PARALLEL_VERBOSE)
            flint_printf("parallel_do with num_threads = %wd\n", num_threads);

        if (num_workers < 1)
        {
            flint_give_back_threads(handles, num_workers);
            return basecase(res, vec, n, ctx);
        }
        else
        {
            work_chunk_t * work;
            slong chunk_size;
            gr_ptr results;

            TMP_INIT;
            TMP_START;

            work = TMP_ALLOC(num_threads * sizeof(work_chunk_t));
            results = TMP_ALLOC(num_threads * ctx->sizeof_elem);

            _gr_vec_init(results, num_threads, ctx);

            chunk_size = (n + num_threads - 1) / num_threads;

            for (i = 0; i < num_threads; i++)
            {
                work[i].f = basecase;
                work[i].vec = vec;
                work[i].res = GR_ENTRY(results, i, ctx->sizeof_elem);
                work[i].a = i * chunk_size;
                work[i].b = FLINT_MIN((i + 1) * chunk_size, n);
                work[i].step = 1;
                work[i].ctx = ctx;
                work[i].status = GR_SUCCESS;
            }

            if ((flags & FLINT_PARALLEL_VERBOSE))
            {
                for (i = 0; i < num_threads; i++)
                {
                    flint_printf("thread #%wd allocated a = %wd, b = %wd, step = %wd\n", i, work[i].a, work[i].b, work[i].step);
                }
            }

            for (i = 0; i < num_workers; i++)
                thread_pool_wake(global_thread_pool, handles[i], 0, worker, &work[i]);

            worker(&work[num_workers]);

            for (i = 0; i < num_workers; i++)
                thread_pool_wait(global_thread_pool, handles[i]);

            for (i = 0; i < num_workers; i++)
                status |= work[i].status;

            status |= basecase(res, results, num_threads, ctx);

            _gr_vec_clear(results, num_threads, ctx);

            flint_give_back_threads(handles, num_workers);
            TMP_END;
        }

        return status;
    }
}

int
_gr_vec_sum_parallel(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    return _gr_vec_parallel_reduce(res, (gr_method_vec_reduce_op) _gr_vec_sum, vec, len, ctx, -1, FLINT_PARALLEL_UNIFORM);
}

typedef struct
{
    gr_srcptr vec;
    gr_ctx_struct * ctx;
}
bsplit_args_t;

typedef struct
{
    gr_ptr res;
    int status;
}
bsplit_res_t;

static void
bsplit_init(bsplit_res_t * x, bsplit_args_t * args)
{
    x->res = gr_heap_init(args->ctx);
    x->status = GR_SUCCESS;
}

static void
bsplit_clear(bsplit_res_t * x, bsplit_args_t * args)
{
    gr_heap_clear(x->res, args->ctx);
}

static void
bsplit_basecase(bsplit_res_t * res, slong a, slong b, bsplit_args_t * args)
{
    res->status = _gr_vec_sum(res->res, GR_ENTRY(args->vec, a, args->ctx->sizeof_elem), b - a, args->ctx);
}

static void
bsplit_merge(bsplit_res_t * res, bsplit_res_t * left, bsplit_res_t * right, bsplit_args_t * args)
{
    res->status = left->status | right->status;

    if (res->status == GR_SUCCESS)
        res->status |= gr_add(res->res, left->res, right->res, args->ctx);
}

int
_gr_vec_sum_bsplit_parallel(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)
{
    bsplit_args_t args;
    bsplit_res_t _res;

    args.vec = vec;
    args.ctx = ctx;

    _res.res = res;
    _res.status = 0;

    flint_parallel_binary_splitting(&_res,
        (bsplit_basecase_func_t) bsplit_basecase,
        (bsplit_merge_func_t) bsplit_merge,
        sizeof(bsplit_res_t),
        (bsplit_init_func_t) bsplit_init,
        (bsplit_clear_func_t) bsplit_clear,
        &args, 0, len, basecase_cutoff, -1, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);

    return _res.status;
}

int
_gr_vec_sum_bsplit(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)
{
    if (len < basecase_cutoff)
    {
        return _gr_vec_sum(res, vec, len, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        status |= _gr_vec_sum_bsplit(res, vec, len / 2, basecase_cutoff, ctx);
        status |= _gr_vec_sum_bsplit(t, GR_ENTRY(vec, len / 2, ctx->sizeof_elem), len - len / 2, basecase_cutoff, ctx);
        status |= gr_add(res, res, t, ctx);

        GR_TMP_CLEAR(t, ctx);
        return status;
    }
}

int
_gr_vec_sum_serial(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;

    if (len <= 2)
    {
        if (len == 2)
            return add(res, vec, GR_ENTRY(vec, 1, ctx->sizeof_elem), ctx);
        else if (len == 1)
            return gr_set(res, vec, ctx);
        else
            return gr_zero(res, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    status |= add(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len; i++)
        status |= add(res, res, GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
_gr_vec_sum_generic(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;

    if (len <= 2)
    {
        if (len == 2)
            return add(res, vec, GR_ENTRY(vec, 1, ctx->sizeof_elem), ctx);
        else if (len == 1)
            return gr_set(res, vec, ctx);
        else
            return gr_zero(res, ctx);
    }

    /* Todo: algorithm selection. Binary splitting is needed to get
       optimal complexity when adding fractions. The cutoff should be lower
       than this when we actually have fractions; on the other hand,
       we don't want to introduce overhead for rings with cheap addition. */
    if (len > 100 && gr_ctx_is_finite(ctx) != T_TRUE)
    {
        return _gr_vec_sum_bsplit(res, vec, len, 100, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    status |= add(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len; i++)
        status |= add(res, res, GR_ENTRY(vec, i, sz), ctx);

    return status;
}
