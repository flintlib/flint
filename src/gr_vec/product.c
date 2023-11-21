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

int _gr_vec_parallel_reduce(gr_ptr res, gr_method_vec_reduce_op basecase, gr_srcptr vec, slong n, gr_ctx_t ctx, int thread_limit, int flags);

int
_gr_vec_product_parallel(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    return _gr_vec_parallel_reduce(res, (gr_method_vec_reduce_op) _gr_vec_product, vec, len, ctx, -1, FLINT_PARALLEL_UNIFORM);
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
    res->status = _gr_vec_product(res->res, GR_ENTRY(args->vec, a, args->ctx->sizeof_elem), b - a, args->ctx);
}

static void
bsplit_merge(bsplit_res_t * res, bsplit_res_t * left, bsplit_res_t * right, bsplit_args_t * args)
{
    res->status = left->status | right->status;

    if (res->status == GR_SUCCESS)
        res->status |= gr_mul(res->res, left->res, right->res, args->ctx);
}

int
_gr_vec_product_bsplit_parallel(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)
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
_gr_vec_product_bsplit(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)
{
    if (len < basecase_cutoff)
    {
        return _gr_vec_product(res, vec, len, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        status |= _gr_vec_product_bsplit(res, vec, len / 2, basecase_cutoff, ctx);
        status |= _gr_vec_product_bsplit(t, GR_ENTRY(vec, len / 2, ctx->sizeof_elem), len - len / 2, basecase_cutoff, ctx);
        status |= gr_mul(res, res, t, ctx);

        GR_TMP_CLEAR(t, ctx);
        return status;
    }
}

int
_gr_vec_product_serial(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    int status;
    slong i, sz;

    if (len <= 2)
    {
        if (len == 2)
            return mul(res, vec, GR_ENTRY(vec, 1, ctx->sizeof_elem), ctx);
        else if (len == 1)
            return gr_set(res, vec, ctx);
        else
            return gr_one(res, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    status |= mul(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len; i++)
        status |= mul(res, res, GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
_gr_vec_product_generic(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    int status;
    slong i, sz;

    if (len <= 2)
    {
        if (len == 2)
            return mul(res, vec, GR_ENTRY(vec, 1, ctx->sizeof_elem), ctx);
        else if (len == 1)
            return gr_set(res, vec, ctx);
        else
            return gr_one(res, ctx);
    }

    /* Todo: algorithm selection. */
    if (len > 20 && gr_ctx_is_finite(ctx) != T_TRUE)
    {
        if (len > 500 && gr_ctx_is_threadsafe(ctx) == T_TRUE)
            return _gr_vec_product_bsplit_parallel(res, vec, len, 500, ctx);
        else
            return _gr_vec_product_bsplit(res, vec, len, 20, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    status |= mul(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len; i++)
        status |= mul(res, res, GR_ENTRY(vec, i, sz), ctx);

    return status;
}
