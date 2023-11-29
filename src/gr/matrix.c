/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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

void
matrix_init(gr_mat_t res, gr_ctx_t ctx)
{
    gr_mat_init(res, MATRIX_CTX(ctx)->nrows, MATRIX_CTX(ctx)->ncols, MATRIX_CTX(ctx)->base_ring);
}

int matrix_ctx_write(gr_stream_t out, gr_ctx_t ctx)
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

truth_t matrix_ctx_is_ring(gr_ctx_t ctx)
{
    truth_t shape_ok;

    shape_ok = (!MATRIX_CTX(ctx)->all_sizes && MATRIX_CTX(ctx)->nrows == MATRIX_CTX(ctx)->ncols) ? T_TRUE : T_FALSE;

    if (shape_ok == T_TRUE && MATRIX_CTX(ctx)->nrows == 0)
        return T_TRUE;

    return truth_and(shape_ok, gr_ctx_is_ring(MATRIX_CTX(ctx)->base_ring));
}

/* todo: public */
truth_t gr_ctx_matrix_is_fixed_size(gr_ctx_t ctx)
{
    return (MATRIX_CTX(ctx)->all_sizes) ? T_FALSE : T_TRUE;
}

truth_t
matrix_ctx_is_threadsafe(gr_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(MATRIX_CTX(ctx)->base_ring);
}

void
matrix_clear(gr_mat_t res, gr_ctx_t ctx)
{
    gr_mat_clear(res, MATRIX_CTX(ctx)->base_ring);
}

void
matrix_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)
{
    gr_mat_swap(mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

void
matrix_set_shallow(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    *res = *mat;
}

int
matrix_write(gr_stream_t out, gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_write(out, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_randtest(gr_mat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        _gr_mat_resize(res, n_randint(state, 7), n_randint(state, 7), MATRIX_CTX(ctx)->base_ring);

    return gr_mat_randtest(res, state, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_equal(mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_set(res, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_si(gr_mat_t res, slong v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_si(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_ui(gr_mat_t res, ulong v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_ui(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_fmpz(gr_mat_t res, const fmpz_t v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_fmpz(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_fmpq(gr_mat_t res, const fmpq_t v, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_set_fmpq(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
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
        slong i, j;
        int status;
        slong sz, xsz;

        if (res->r != xmat->r || res->c != xmat->c)
        {
            if (MATRIX_CTX(ctx)->all_sizes)
                _gr_mat_resize(res, xmat->r, xmat->c, MATRIX_CTX(ctx)->base_ring);
            else
                return GR_DOMAIN;
        }

        sz = MATRIX_CTX(ctx)->base_ring->sizeof_elem;
        xsz = MATRIX_CTX(x_ctx)->base_ring->sizeof_elem;

        for (i = 0; i < xmat->r; i++)
        {
            for (j = 0; j < xmat->c; j++)
            {
                status = gr_set_other(GR_MAT_ENTRY(res, i, j, sz),
                            GR_MAT_ENTRY(xmat, i, j, xsz),
                            MATRIX_CTX(x_ctx)->base_ring,
                            MATRIX_CTX(ctx)->base_ring);

                if (status != GR_SUCCESS)
                    return status;
            }
        }

        return GR_SUCCESS;
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

int
matrix_zero(gr_mat_t res, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_zero(res, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_one(gr_mat_t res, gr_ctx_t ctx)
{
    if (MATRIX_CTX(ctx)->all_sizes)
        return GR_DOMAIN;

    return gr_mat_one(res, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_is_zero(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_zero(mat, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_is_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_one(mat, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_neg_one(mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_neg(res, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    if (mat1->r != mat2->r || mat1->c != mat2->c)
        return GR_DOMAIN;

    if (res->r != mat1->r || res->c != mat1->c)
        _gr_mat_resize(res, mat1->r, mat1->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_add(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    if (mat1->r != mat2->r || mat1->c != mat2->c)
        return GR_DOMAIN;

    if (res->r != mat1->r || res->c != mat1->c)
        _gr_mat_resize(res, mat1->r, mat1->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_sub(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
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

    return gr_mat_mul_classical(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

/* todo: fmpz etc */
int
matrix_mul_other(gr_mat_t res, const gr_mat_t mat, gr_ptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    if (y_ctx == ctx)
    {
        return matrix_mul(res, mat, y, ctx);
    }
    else if (y_ctx == MATRIX_CTX(ctx)->base_ring)
    {
        int status = GR_SUCCESS;

        if (res->r != mat->r || res->c != mat->c)
            status = _gr_mat_check_resize(res, mat->r, mat->c, ctx);

        if (status != GR_SUCCESS)
            return status;

        /* todo: resizing */
        return gr_mat_mul_scalar(res, mat, y, y_ctx);
    }
    else if (y_ctx->which_ring == GR_CTX_GR_MAT)
    {
        int status;
        gr_mat_t tmp;
        gr_mat_init(tmp, ((gr_mat_struct *) y)->r, ((gr_mat_struct *) y)->c, MATRIX_CTX(ctx)->base_ring);
        status = matrix_set_other(tmp, y, y_ctx, ctx);
        if (status == GR_SUCCESS)
        {
            status = matrix_mul(res, mat, tmp, ctx);
        }
        gr_mat_clear(tmp, MATRIX_CTX(ctx)->base_ring);
        return status;
    }
    else
    {
        int status;
        gr_ptr c;

        GR_TMP_INIT(c, MATRIX_CTX(ctx)->base_ring);

        status = gr_set_other(c, y, y_ctx, MATRIX_CTX(ctx)->base_ring);

        if (status == GR_SUCCESS)
        {
            if (res->r != mat->r || res->c != mat->c)
                status = _gr_mat_check_resize(res, mat->r, mat->c, ctx);

            if (status == GR_SUCCESS)
                status = gr_mat_mul_scalar(res, mat, c, MATRIX_CTX(ctx)->base_ring);
        }

        GR_TMP_CLEAR(c, MATRIX_CTX(ctx)->base_ring);

        return status;
    }

    return GR_UNABLE;
}

/* todo: fmpz etc */
int
matrix_div_other(gr_mat_t res, const gr_mat_t mat, gr_ptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    /* todo: matrix-matrix division ? */

    if (y_ctx == MATRIX_CTX(ctx)->base_ring)
    {
        int status = GR_SUCCESS;

        if (res->r != mat->r || res->c != mat->c)
            status = _gr_mat_check_resize(res, mat->r, mat->c, ctx);

        if (status != GR_SUCCESS)
            return status;

        return gr_mat_div_scalar(res, mat, y, y_ctx);
    }
    else
    {
        int status;
        gr_ptr c;

        GR_TMP_INIT(c, MATRIX_CTX(ctx)->base_ring);

        status = gr_set_other(c, y, y_ctx, MATRIX_CTX(ctx)->base_ring);

        if (status == GR_SUCCESS)
        {
            if (res->r != mat->r || res->c != mat->c)
                status = _gr_mat_check_resize(res, mat->r, mat->c, ctx);

            if (status == GR_SUCCESS)
                status = gr_mat_div_scalar(res, mat, c, MATRIX_CTX(ctx)->base_ring);
        }

        GR_TMP_CLEAR(c, MATRIX_CTX(ctx)->base_ring);

        return status;
    }

    return GR_UNABLE;
}

int
matrix_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (mat->r != mat->c)
        return GR_DOMAIN;

    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->r, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_inv(res, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_exp(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_exp(res, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_log(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    if (res->r != mat->r || res->c != mat->c)
        _gr_mat_resize(res, mat->r, mat->c, MATRIX_CTX(ctx)->base_ring);

    return gr_mat_log(res, mat, MATRIX_CTX(ctx)->base_ring);
}


int _gr_mat_methods_initialized = 0;

gr_static_method_table _gr_mat_methods;

gr_method_tab_input _gr_mat_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) matrix_ctx_write},
    {GR_METHOD_CTX_IS_RING, (gr_funcptr) matrix_ctx_is_ring},
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
    {GR_METHOD_SUB,         (gr_funcptr) matrix_sub},
    {GR_METHOD_MUL,         (gr_funcptr) matrix_mul},
    {GR_METHOD_MUL_OTHER,   (gr_funcptr) matrix_mul_other},
    {GR_METHOD_DIV_OTHER,   (gr_funcptr) matrix_div_other},
    {GR_METHOD_INV,         (gr_funcptr) matrix_inv},
    {GR_METHOD_EXP,         (gr_funcptr) matrix_exp},
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
