/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Matrices over generic rings */

#include "fmpz.h"
#include "fmpq.h"
#include "gr.h"
#include "gr_mat.h"

/* todo: for matrix "ring", verify that element domain is a ring */

/* todo: recycle storage? */
void
_gr_mat_resize(gr_mat_t mat, slong r, slong c, gr_ctx_t ctx)
{
    gr_mat_clear(mat, ctx);
    gr_mat_init(mat, r, c, ctx);
}

int
_gr_mat_check_resize(gr_mat_t mat, slong r, slong c, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
    {
        _gr_mat_resize(mat, r, c, MATRIX_CTX(ctx)->base_ring);
        return GR_SUCCESS;
    }
    else
    {
        if (r != MATRIX_CTX(ctx)->nrows || c != MATRIX_CTX(ctx)->ncols)
            return GR_DOMAIN;

        if (mat->r != r || mat->c != c)
            _gr_mat_resize(mat, r, c, MATRIX_CTX(ctx)->base_ring);

        return GR_SUCCESS;
    }
}

static void
matrix_init(gr_mat_t res, gr_ctx_t ctx)
{
    gr_mat_init(res, MATRIX_CTX(ctx)->nrows, MATRIX_CTX(ctx)->ncols, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_ctx_ptr elem_ctx = MATRIX_CTX(ctx)->base_ring;

    if (MATRIX_CTX(ctx)->all_sizes)
    {
        gr_stream_write(out, "Matrices (any shape) over ");
    }
    else
    {
        if (gr_ctx_is_ring(ctx) == T_TRUE)
            gr_stream_write(out, "Ring of ");
        else
            gr_stream_write(out, "Space of ");

        gr_stream_write_si(out, MATRIX_CTX(ctx)->nrows);
        gr_stream_write(out, " x ");
        gr_stream_write_si(out, MATRIX_CTX(ctx)->ncols);
        gr_stream_write(out, " ");

        gr_stream_write(out, "matrices over ");
    }

    gr_ctx_write(out, elem_ctx);
    return GR_SUCCESS;
}

static truth_t matrix_ctx_is_ring(gr_ctx_t ctx)
{
    int shape_ok = (!MATRIX_CTX(ctx)->all_sizes && MATRIX_CTX(ctx)->nrows == MATRIX_CTX(ctx)->ncols);

    if (!shape_ok)
        return T_FALSE;

    if (MATRIX_CTX(ctx)->nrows == 0)
        return T_TRUE;

    return gr_ctx_is_ring(MATRIX_CTX(ctx)->base_ring);
}

static truth_t matrix_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    int shape_ok = (!MATRIX_CTX(ctx)->all_sizes && MATRIX_CTX(ctx)->nrows == MATRIX_CTX(ctx)->ncols);

    if (!shape_ok)
        return T_FALSE;

    if (MATRIX_CTX(ctx)->nrows == 0)
        return T_TRUE;

    if (MATRIX_CTX(ctx)->nrows == 1)
        return gr_ctx_is_commutative_ring(MATRIX_CTX(ctx)->base_ring);

    return gr_ctx_is_zero_ring(MATRIX_CTX(ctx)->base_ring);
}

/* todo: public */
truth_t gr_ctx_matrix_is_fixed_size(gr_ctx_t ctx)
{
    return (MATRIX_CTX(ctx)->all_sizes) ? T_FALSE : T_TRUE;
}

static truth_t
matrix_ctx_is_threadsafe(gr_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(MATRIX_CTX(ctx)->base_ring);
}

static void
matrix_clear(gr_mat_t res, gr_ctx_t ctx)
{
    gr_mat_clear(res, MATRIX_CTX(ctx)->base_ring);
}

static void
matrix_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)
{
    gr_mat_swap(mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

static void
matrix_set_shallow(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    *res = *mat;
}

static int
matrix_write(gr_stream_t out, gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_write(out, mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_randtest(gr_mat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        _gr_mat_resize(res, n_randint(state, 7), n_randint(state, 7), MATRIX_CTX(ctx)->base_ring);

    return gr_mat_randtest(res, state, MATRIX_CTX(ctx)->base_ring);
}

static truth_t
matrix_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_equal(mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_set(res, mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_set_si(gr_mat_t res, slong v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_si(res, v, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_set_ui(gr_mat_t res, ulong v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_ui(res, v, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_set_fmpz(gr_mat_t res, const fmpz_t v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_fmpz(res, v, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_set_fmpq(gr_mat_t res, const fmpq_t v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_fmpq(res, v, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_set_other(gr_mat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return matrix_set(res, x, ctx);
    }
    else if (x_ctx == MATRIX_CTX(ctx)->base_ring)
    {
        if (MATRIX_CTX(ctx)->all_sizes)
            return GR_DOMAIN;

        return gr_mat_set_scalar(res, x, x_ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_MAT)
    {
        const gr_mat_struct * xmat = x;

        if (res->r != xmat->r || res->c != xmat->c)
        {
            if (MATRIX_CTX(ctx)->all_sizes)
                _gr_mat_resize(res, xmat->r, xmat->c, MATRIX_CTX(ctx)->base_ring);
            else
                return GR_DOMAIN;
        }

        return gr_mat_set_gr_mat_other(res, xmat, MATRIX_CTX(x_ctx)->base_ring,
                            MATRIX_CTX(ctx)->base_ring);
    }
    else
    {
        int status = GR_SUCCESS;
        gr_ptr tmp;

        if (MATRIX_CTX(ctx)->all_sizes)
            return GR_UNABLE;

        GR_TMP_INIT(tmp, MATRIX_CTX(ctx)->base_ring);

        status = gr_set_other(tmp, x, x_ctx, MATRIX_CTX(ctx)->base_ring);

        if (status == GR_SUCCESS)
            status = gr_mat_set_scalar(res, tmp, MATRIX_CTX(ctx)->base_ring);

        GR_TMP_CLEAR(tmp, MATRIX_CTX(ctx)->base_ring);
        return status;
    }
}

static int
matrix_zero(gr_mat_t res, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_zero(res, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_one(gr_mat_t res, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_one(res, MATRIX_CTX(ctx)->base_ring);
}

static truth_t
matrix_is_zero(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_zero(mat, MATRIX_CTX(ctx)->base_ring);
}

static truth_t
matrix_is_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_one(mat, MATRIX_CTX(ctx)->base_ring);
}

static truth_t
matrix_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_neg_one(mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_neg(res, mat, MATRIX_CTX(ctx)->base_ring);
}

#define DEF_MAT_OP_1(op, gr_mat_op) static int \
    op (gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx) \
    { \
        if (mat1->r != mat2->r || mat1->c != mat2->c) \
            return GR_DOMAIN; \
        if (res->r != mat1->r || res->c != mat1->c) \
            _gr_mat_resize(res, mat1->r, mat1->c, MATRIX_CTX(ctx)->base_ring); \
        return gr_mat_op(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring); \
    }

#define DEF_MAT_SCALAR_OP_1(op, typ, gr_mat_op_scalar) static int \
    op (gr_mat_t res, const gr_mat_t mat1, typ x, gr_ctx_t ctx) \
    { \
        if (res->r != mat1->r || res->c != mat1->c) \
            _gr_mat_resize(res, mat1->r, mat1->c, MATRIX_CTX(ctx)->base_ring); \
        return gr_mat_op_scalar(res, mat1, x, MATRIX_CTX(ctx)->base_ring); \
    }

#define DEF_MAT_OP_OTHER(op, matrix_op, gr_mat_op_scalar) static int \
    op(gr_mat_t res, const gr_mat_t mat, gr_ptr y, gr_ctx_t y_ctx, gr_ctx_t ctx) \
    { \
        if (y_ctx == ctx) \
        { \
            return matrix_op(res, mat, y, ctx); \
        } \
        else if (y_ctx == MATRIX_CTX(ctx)->base_ring) \
        { \
            int status = GR_SUCCESS; \
            if (res->r != mat->r || res->c != mat->c) \
                status = _gr_mat_check_resize(res, mat->r, mat->c, ctx); \
            if (status != GR_SUCCESS) \
                return status; \
            return gr_mat_op_scalar(res, mat, y, y_ctx); \
        } \
        else if (y_ctx->which_ring == GR_CTX_GR_MAT) \
        { \
            int status; \
            gr_mat_t tmp; \
            gr_mat_init(tmp, ((gr_mat_struct *) y)->r, ((gr_mat_struct *) y)->c, MATRIX_CTX(ctx)->base_ring); \
            status = matrix_set_other(tmp, y, y_ctx, ctx); \
            if (status == GR_SUCCESS) \
            { \
                status = matrix_op(res, mat, tmp, ctx); \
            } \
            gr_mat_clear(tmp, MATRIX_CTX(ctx)->base_ring); \
            return status; \
        } \
        else \
        { \
            int status; \
            gr_ptr c; \
            GR_TMP_INIT(c, MATRIX_CTX(ctx)->base_ring); \
            status = gr_set_other(c, y, y_ctx, MATRIX_CTX(ctx)->base_ring); \
            if (status == GR_SUCCESS) \
            { \
                if (res->r != mat->r || res->c != mat->c) \
                    status = _gr_mat_check_resize(res, mat->r, mat->c, ctx); \
                if (status == GR_SUCCESS) \
                    status = gr_mat_op_scalar(res, mat, c, MATRIX_CTX(ctx)->base_ring); \
            } \
            GR_TMP_CLEAR(c, MATRIX_CTX(ctx)->base_ring); \
            return status; \
        } \
        return GR_UNABLE; \
    } \

#define DEF_MAT_OTHER_OP(op, matrix_op, gr_scalar_op_mat) static int \
    op(gr_mat_t res, gr_ptr x, gr_ctx_t x_ctx, const gr_mat_t mat, gr_ctx_t ctx) \
    { \
        if (x_ctx == ctx) \
        { \
            return matrix_op(res, x, mat, ctx); \
        } \
        else if (x_ctx == MATRIX_CTX(ctx)->base_ring) \
        { \
            int status = GR_SUCCESS; \
            if (res->r != mat->r || res->c != mat->c) \
                status = _gr_mat_check_resize(res, mat->r, mat->c, ctx); \
            if (status != GR_SUCCESS) \
                return status; \
            return gr_scalar_op_mat(res, x, mat, x_ctx); \
        } \
        else if (x_ctx->which_ring == GR_CTX_GR_MAT) \
        { \
            int status; \
            gr_mat_t tmp; \
            gr_mat_init(tmp, ((gr_mat_struct *) x)->r, ((gr_mat_struct *) x)->c, MATRIX_CTX(ctx)->base_ring); \
            status = matrix_set_other(tmp, x, x_ctx, ctx); \
            if (status == GR_SUCCESS) \
            { \
                status = matrix_op(res, tmp, mat, ctx); \
            } \
            gr_mat_clear(tmp, MATRIX_CTX(ctx)->base_ring); \
            return status; \
        } \
        else \
        { \
            int status; \
            gr_ptr c; \
            GR_TMP_INIT(c, MATRIX_CTX(ctx)->base_ring); \
            status = gr_set_other(c, x, x_ctx, MATRIX_CTX(ctx)->base_ring); \
            if (status == GR_SUCCESS) \
            { \
                if (res->r != mat->r || res->c != mat->c) \
                    status = _gr_mat_check_resize(res, mat->r, mat->c, ctx); \
                if (status == GR_SUCCESS) \
                    status = gr_scalar_op_mat(res, c, mat, MATRIX_CTX(ctx)->base_ring); \
            } \
            GR_TMP_CLEAR(c, MATRIX_CTX(ctx)->base_ring); \
            return status; \
        } \
        return GR_UNABLE; \
    } \

DEF_MAT_OP_1(matrix_add, gr_mat_add)
DEF_MAT_OP_1(matrix_sub, gr_mat_sub)

static int
matrix_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    if (mat1->c != mat2->r)
        return GR_DOMAIN;

    if (res->r != mat1->r || res->c != mat2->c)
    {
        if (res == mat1 || res == mat2)
        {
            int status;
            gr_mat_t tmp;
            gr_mat_init(tmp, mat1->r, mat2->c, MATRIX_CTX(ctx)->base_ring);
            status = matrix_mul(tmp, mat1, mat2, ctx);
            gr_mat_swap(res, tmp, MATRIX_CTX(ctx)->base_ring);
            gr_mat_clear(tmp, MATRIX_CTX(ctx)->base_ring);
            return status;
        }
        else
        {
            _gr_mat_resize(res, mat1->r, mat2->c, MATRIX_CTX(ctx)->base_ring);
        }
    }

    return gr_mat_mul(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_div(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

static int
matrix_pow(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

DEF_MAT_OP_OTHER(matrix_add_other, matrix_add, gr_mat_add_scalar)
DEF_MAT_OP_OTHER(matrix_sub_other, matrix_sub, gr_mat_sub_scalar)
DEF_MAT_OP_OTHER(matrix_mul_other, matrix_mul, gr_mat_mul_scalar)
DEF_MAT_OP_OTHER(matrix_div_other, matrix_div, gr_mat_div_scalar)
DEF_MAT_OP_OTHER(matrix_pow_other, matrix_pow, gr_mat_pow_scalar)

DEF_MAT_OTHER_OP(matrix_other_add, matrix_add, gr_mat_scalar_add)
DEF_MAT_OTHER_OP(matrix_other_sub, matrix_sub, gr_mat_scalar_sub)
DEF_MAT_OTHER_OP(matrix_other_mul, matrix_mul, gr_mat_scalar_mul)

DEF_MAT_SCALAR_OP_1(matrix_add_ui, ulong, gr_mat_add_ui)
DEF_MAT_SCALAR_OP_1(matrix_add_si, slong, gr_mat_add_si)
DEF_MAT_SCALAR_OP_1(matrix_add_fmpz, const fmpz_t, gr_mat_add_fmpz)
DEF_MAT_SCALAR_OP_1(matrix_add_fmpq, const fmpq_t, gr_mat_add_fmpq)

DEF_MAT_SCALAR_OP_1(matrix_sub_ui, ulong, gr_mat_sub_ui)
DEF_MAT_SCALAR_OP_1(matrix_sub_si, slong, gr_mat_sub_si)
DEF_MAT_SCALAR_OP_1(matrix_sub_fmpz, const fmpz_t, gr_mat_sub_fmpz)
DEF_MAT_SCALAR_OP_1(matrix_sub_fmpq, const fmpq_t, gr_mat_sub_fmpq)

DEF_MAT_SCALAR_OP_1(matrix_mul_ui, ulong, gr_mat_mul_ui)
DEF_MAT_SCALAR_OP_1(matrix_mul_si, slong, gr_mat_mul_si)
DEF_MAT_SCALAR_OP_1(matrix_mul_fmpz, const fmpz_t, gr_mat_mul_fmpz)
DEF_MAT_SCALAR_OP_1(matrix_mul_fmpq, const fmpq_t, gr_mat_mul_fmpq)

DEF_MAT_SCALAR_OP_1(matrix_div_ui, ulong, gr_mat_div_ui)
DEF_MAT_SCALAR_OP_1(matrix_div_si, slong, gr_mat_div_si)
DEF_MAT_SCALAR_OP_1(matrix_div_fmpz, const fmpz_t, gr_mat_div_fmpz)
DEF_MAT_SCALAR_OP_1(matrix_div_fmpq, const fmpq_t, gr_mat_div_fmpq)

DEF_MAT_SCALAR_OP_1(matrix_pow_ui, ulong, gr_mat_pow_ui)
DEF_MAT_SCALAR_OP_1(matrix_pow_si, slong, gr_mat_pow_si)
DEF_MAT_SCALAR_OP_1(matrix_pow_fmpz, const fmpz_t, gr_mat_pow_fmpz)
DEF_MAT_SCALAR_OP_1(matrix_pow_fmpq, const fmpq_t, gr_mat_pow_fmpq)


static int
matrix_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (mat->r != mat->c)
        return GR_DOMAIN;

    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->r, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_inv(res, mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_exp(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_exp(res, mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_log(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_log(res, mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_sqrt(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_sqrt(res, mat, MATRIX_CTX(ctx)->base_ring);
}

static int
matrix_rsqrt(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_rsqrt(res, mat, MATRIX_CTX(ctx)->base_ring);
}


int _gr_mat_methods_initialized = 0;

gr_static_method_table _gr_mat_methods;

gr_method_tab_input _gr_mat_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) matrix_ctx_write},
    {GR_METHOD_CTX_IS_RING, (gr_funcptr) matrix_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) matrix_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) matrix_ctx_is_threadsafe},
    {GR_METHOD_INIT,        (gr_funcptr) matrix_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) matrix_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) matrix_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) matrix_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) matrix_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) matrix_write},
    {GR_METHOD_ZERO,        (gr_funcptr) matrix_zero},
    {GR_METHOD_ONE,         (gr_funcptr) matrix_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) matrix_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) matrix_is_one},
    {GR_METHOD_IS_NEG_ONE,  (gr_funcptr) matrix_is_neg_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) matrix_equal},
    {GR_METHOD_SET,         (gr_funcptr) matrix_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) matrix_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) matrix_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) matrix_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) matrix_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) matrix_set_other},
    {GR_METHOD_NEG,         (gr_funcptr) matrix_neg},
    {GR_METHOD_ADD,         (gr_funcptr) matrix_add},
    {GR_METHOD_ADD_OTHER,   (gr_funcptr) matrix_add_other},
    {GR_METHOD_OTHER_ADD,   (gr_funcptr) matrix_other_add},
    {GR_METHOD_ADD_UI,      (gr_funcptr) matrix_add_ui},
    {GR_METHOD_ADD_SI,      (gr_funcptr) matrix_add_si},
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) matrix_add_fmpz},
    {GR_METHOD_ADD_FMPQ,    (gr_funcptr) matrix_add_fmpq},
    {GR_METHOD_SUB,         (gr_funcptr) matrix_sub},
    {GR_METHOD_SUB_OTHER,   (gr_funcptr) matrix_sub_other},
    {GR_METHOD_OTHER_SUB,   (gr_funcptr) matrix_other_sub},
    {GR_METHOD_SUB_UI,      (gr_funcptr) matrix_sub_ui},
    {GR_METHOD_SUB_SI,      (gr_funcptr) matrix_sub_si},
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) matrix_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,    (gr_funcptr) matrix_sub_fmpq},
    {GR_METHOD_MUL,         (gr_funcptr) matrix_mul},
    {GR_METHOD_MUL_OTHER,   (gr_funcptr) matrix_mul_other},
    {GR_METHOD_OTHER_MUL,   (gr_funcptr) matrix_other_mul},
    {GR_METHOD_MUL_UI,      (gr_funcptr) matrix_mul_ui},
    {GR_METHOD_MUL_SI,      (gr_funcptr) matrix_mul_si},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) matrix_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) matrix_mul_fmpq},
    {GR_METHOD_DIV_OTHER,   (gr_funcptr) matrix_div_other},
    {GR_METHOD_DIV_UI,      (gr_funcptr) matrix_div_ui},
    {GR_METHOD_DIV_SI,      (gr_funcptr) matrix_div_si},
    {GR_METHOD_DIV_FMPZ,    (gr_funcptr) matrix_div_fmpz},
    {GR_METHOD_DIV_FMPQ,    (gr_funcptr) matrix_div_fmpq},
    {GR_METHOD_POW_OTHER,   (gr_funcptr) matrix_pow_other},
    {GR_METHOD_POW_UI,      (gr_funcptr) matrix_pow_ui},
    {GR_METHOD_POW_SI,      (gr_funcptr) matrix_pow_si},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) matrix_pow_fmpz},
    {GR_METHOD_POW_FMPQ,    (gr_funcptr) matrix_pow_fmpq},
    {GR_METHOD_INV,         (gr_funcptr) matrix_inv},
    {GR_METHOD_EXP,         (gr_funcptr) matrix_exp},
    {GR_METHOD_SQRT,        (gr_funcptr) matrix_sqrt},
    {GR_METHOD_RSQRT,       (gr_funcptr) matrix_rsqrt},
    {GR_METHOD_LOG,         (gr_funcptr) matrix_log},
    {0,                     (gr_funcptr) NULL},
};

void
_gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, int all_sizes, slong nrows, slong ncols)
{
    ctx->which_ring = GR_CTX_GR_MAT;
    ctx->sizeof_elem = sizeof(gr_mat_struct);
    ctx->size_limit = WORD_MAX;

    if (nrows < 0 || ncols < 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    MATRIX_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    MATRIX_CTX(ctx)->all_sizes = all_sizes;
    MATRIX_CTX(ctx)->nrows = nrows;
    MATRIX_CTX(ctx)->ncols = ncols;

    ctx->methods = _gr_mat_methods;

    if (!_gr_mat_methods_initialized)
    {
        gr_method_tab_init(_gr_mat_methods, _gr_mat_methods_input);
        _gr_mat_methods_initialized = 1;
    }
}

void gr_ctx_init_matrix_domain(gr_ctx_t ctx, gr_ctx_t base_ring)
{
    _gr_ctx_init_matrix(ctx, base_ring, 1, 0, 0);
}

void gr_ctx_init_matrix_space(gr_ctx_t ctx, gr_ctx_t base_ring, slong nrows, slong ncols)
{
    _gr_ctx_init_matrix(ctx, base_ring, 0, nrows, ncols);
}
