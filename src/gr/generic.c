/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arith.h"
#include "bernoulli.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "gr_special.h"

#define DEBUG_RINGS 0

/* todo: rings only */
#if DEBUG_RINGS
void gr_generic_init(void) { flint_printf("ctx must implement init()\n"); flint_abort(); }
void gr_generic_clear(void) { flint_printf("ctx must implement clear()\n"); flint_abort(); }
void gr_generic_swap(void) { flint_printf("ctx must implement swap()\n"); flint_abort(); }
void gr_generic_randtest(void) { flint_printf("ctx must implement randtest()\n"); flint_abort(); }
void gr_generic_write(void) { flint_printf("ctx must implement write()\n"); flint_abort(); }
void gr_generic_zero(void) { flint_printf("ctx must implement zero()\n"); flint_abort(); }
void gr_generic_one(void) { flint_printf("ctx must implement one()\n"); flint_abort(); }
void gr_generic_equal(void) { flint_printf("ctx must implement equal()\n"); flint_abort(); }
void gr_generic_set(void) { flint_printf("ctx must implement set()\n"); flint_abort(); }
void gr_generic_set_si(void) { flint_printf("ctx must implement set_si()\n"); flint_abort(); }
void gr_generic_set_ui(void) { flint_printf("ctx must implement set_ui()\n"); flint_abort(); }
void gr_generic_set_fmpz(void) { flint_printf("ctx must implement set_fmpz()\n"); flint_abort(); }
void gr_generic_neg(void) { flint_printf("ctx must implement neg()\n"); flint_abort(); }
void gr_generic_add(void) { flint_printf("ctx must implement add()\n"); flint_abort(); }
void gr_generic_sub(void) { flint_printf("ctx must implement sub()\n"); flint_abort(); }
void gr_generic_mul(void) { flint_printf("ctx must implement mul()\n"); flint_abort(); }
#else
#define gr_generic_init gr_not_implemented
#define gr_generic_clear gr_not_implemented
#define gr_generic_swap gr_not_implemented
#define gr_generic_randtest gr_not_implemented
#define gr_generic_write gr_not_implemented
#define gr_generic_zero gr_not_implemented
#define gr_generic_one gr_not_implemented
#define gr_generic_equal gr_not_implemented
#define gr_generic_set gr_not_implemented
#define gr_generic_set_si gr_not_implemented
#define gr_generic_set_ui gr_not_implemented
#define gr_generic_set_fmpz gr_not_implemented
#define gr_generic_neg gr_not_implemented
#define gr_generic_add gr_not_implemented
#define gr_generic_sub gr_not_implemented
#define gr_generic_mul gr_not_implemented
#endif

truth_t gr_generic_ctx_predicate(gr_ctx_t ctx)
{
    return T_UNKNOWN;
}

truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx)
{
    return T_TRUE;
}

truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx)
{
    return T_FALSE;
}

void
gr_generic_set_shallow(gr_ptr res, gr_srcptr x, const gr_ctx_t ctx)
{
    memcpy(res, x, ctx->sizeof_elem);
}


int gr_generic_randtest_not_zero(gr_ptr x, flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    truth_t is_zero;
    int status = GR_SUCCESS;

    for (i = 0; i < 5; i++)
    {
        status |= gr_randtest(x, state, ctx);

        is_zero = gr_is_zero(x, ctx);
        if (is_zero == T_FALSE)
            return GR_SUCCESS;
    }

    if (n_randint(state, 2))
        status |= gr_one(x, ctx);
    else
        status |= gr_neg_one(x, ctx);

    /* unused */
    status = status;

    is_zero = gr_is_zero(x, ctx);
    if (is_zero == T_FALSE)
        return GR_SUCCESS;

    /* We are in the zero ring */
    if (is_zero == T_TRUE)
        return GR_DOMAIN;

    return GR_UNABLE;
}

/* Generic arithmetic functions */

truth_t gr_generic_is_zero(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t eq;

    GR_TMP_INIT(t, ctx);
    eq = gr_equal(x, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    return eq;
}

truth_t gr_generic_is_one(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t eq;

    GR_TMP_INIT(t, ctx);

    if (gr_one(t, ctx) == GR_SUCCESS)
        eq = gr_equal(x, t, ctx);
    else
        eq = T_UNKNOWN;

    GR_TMP_CLEAR(t, ctx);

    return eq;
}

truth_t gr_generic_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t eq;

    GR_TMP_INIT(t, ctx);

    if (gr_neg_one(t, ctx) == GR_SUCCESS)
        eq = gr_equal(x, t, ctx);
    else
        eq = T_UNKNOWN;

    GR_TMP_CLEAR(t, ctx);

    return eq;
}

int gr_generic_neg_one(gr_ptr res, gr_ctx_t ctx)
{
    int status;
    status = gr_one(res, ctx);
    status |= gr_neg(res, res, ctx);
    return status;
}

int gr_generic_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t xctx, gr_ctx_t ctx)
{
    if (xctx == ctx)
    {
        return gr_set(res, x, ctx);
    }
    else if (xctx->which_ring == GR_CTX_FMPZ)
    {
        return gr_set_fmpz(res, x, ctx);
    }
    else if (xctx->which_ring == GR_CTX_FMPQ)
    {
        return gr_set_fmpq(res, x, ctx);
    }
    else
    {
        return GR_UNABLE;
    }
}

int gr_generic_set_fmpq(gr_ptr res, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t, u;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_set_fmpz(t, fmpq_numref(y), ctx);
    status |= gr_set_fmpz(u, fmpq_denref(y), ctx);

    if (status == GR_SUCCESS)
        status = gr_inv(u, u, ctx);

    if (status == GR_SUCCESS)
        status = gr_mul(res, t, u, ctx);

    GR_TMP_CLEAR2(t, u, ctx);
    return status;
}

int gr_generic_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_fmpz(t, y, ctx);

    if (status == GR_SUCCESS)
        status = gr_add(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_ui(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_si(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_fmpq(t, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_add(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, y, y_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_add(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, x, x_ctx, ctx);

    if (status == GR_SUCCESS)
        status = gr_add(res, t, y, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_neg_ui(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, y);
    fmpz_neg(t, t);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    int status;
    fmpz_t t;
    fmpz_init(t);
    fmpz_neg(t, y);
    status = gr_add_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    int status;
    fmpq_t t;
    fmpq_init(t);
    fmpq_neg(t, y);
    status = gr_add_fmpq(res, x, t, ctx);
    fmpq_clear(t);
    return status;
}

int gr_generic_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, y, y_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_sub(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, x, x_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_sub(res, t, y, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_fmpz(t, y, ctx);

    if (status == GR_SUCCESS)
        status = gr_mul(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_ui(t, y);
    status = gr_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_si(t, y);
    status = gr_mul_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_fmpq(t, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, y, y_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, x, x_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_mul(res, t, y, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul(t, x, y, ctx);
    status |= gr_add(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_ui(t, x, y, ctx);
    status |= gr_add(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_si(t, x, y, ctx);
    status |= gr_add(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_fmpz(t, x, y, ctx);
    status |= gr_add(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_fmpq(t, x, y, ctx);
    status |= gr_add(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, y, y_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_addmul(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul(t, x, y, ctx);
    status |= gr_sub(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_ui(t, x, y, ctx);
    status |= gr_sub(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_si(t, x, y, ctx);
    status |= gr_sub(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_fmpz(t, x, y, ctx);
    status |= gr_sub(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_mul_fmpq(t, x, y, ctx);
    status |= gr_sub(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, y, y_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_submul(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}


int gr_generic_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return gr_add(res, x, x, ctx);
}

int gr_generic_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return gr_mul(res, x, x, ctx);
}

int gr_generic_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    if (y == 0)
    {
        return gr_set(res, x, ctx);
    }
    else
    {
        gr_ptr t;
        int status;

        GR_TMP_INIT(t, ctx);

        status = gr_set_ui(t, 2, ctx);

        if (y >= 0)
        {
            status |= gr_pow_ui(t, t, y, ctx);
            status |= gr_mul(res, x, t, ctx);
        }
        else
        {
            status |= gr_pow_ui(t, t, -y, ctx);
            status |= gr_div(res, x, t, ctx);
        }

        GR_TMP_CLEAR(t, ctx);

        return status;
    }
}

int gr_generic_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return gr_set(res, x, ctx);
    }
    else
    {
        gr_ptr t;
        int status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        status = gr_set_ui(t, 2, ctx);

        if (fmpz_sgn(y) > 0)
        {
            status |= gr_pow_fmpz(t, t, y, ctx);
            status |= gr_mul(res, x, t, ctx);
        }
        else
        {
            fmpz_t u;
            fmpz_init(u);
            fmpz_neg(u, y);
            status |= gr_pow_fmpz(t, t, u, ctx);
            status |= gr_div(res, x, t, ctx);
            fmpz_clear(u);
        }

        GR_TMP_CLEAR(t, ctx);

        return status;
    }
}

int gr_generic_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return gr_set_fmpz(res, x, ctx);
    }
    else
    {
        int status;

        status = gr_set_ui(res, 2, ctx);
        status |= gr_pow_fmpz(res, res, y, ctx);
        status |= gr_mul_fmpz(res, res, x, ctx);

        return status;
    }
}

int gr_generic_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_ptr x, gr_ctx_t ctx)
{
    int status;

    fmpq_t v;
    fmpq_init(v);

    status = gr_get_fmpq(v, x, ctx);

    if (status == GR_SUCCESS)
    {
        slong nbits, dbits;

        dbits = fmpz_val2(fmpq_denref(v));
        fmpz_tdiv_q_2exp(fmpq_denref(v), fmpq_denref(v), dbits);

        if (fmpz_is_one(fmpq_denref(v)))
        {
            nbits = fmpz_val2(fmpq_numref(v));
            fmpz_tdiv_q_2exp(fmpq_numref(v), fmpq_numref(v), nbits);

            fmpz_swap(res1, fmpq_numref(v));
            fmpz_set_si(res2, nbits - dbits);
        }
        else
        {
            status = GR_DOMAIN;
        }
    }

    fmpq_clear(v);

    return status;
}

int gr_generic_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    if (gr_is_neg_one(x, ctx) == T_TRUE)
        return gr_neg_one(res, ctx);

    /* todo: dubious in the zero ring, if comparing with 1 above
       somehow failed */
    if (gr_is_zero(x, ctx) == T_TRUE)
        return GR_DOMAIN;

    return GR_UNABLE;
}

truth_t gr_generic_is_invertible(gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    if (gr_is_neg_one(x, ctx) == T_TRUE)
        return T_TRUE;

    /* todo: dubious in the zero ring, if comparing with 1 above
       somehow failed */
    if (gr_is_zero(x, ctx) == T_TRUE)
        return T_FALSE;

    return T_UNKNOWN;
}


int gr_generic_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_fmpz(t, y, ctx);

    if (status == GR_SUCCESS)
        status = gr_div(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_ui(t, y);
    status = gr_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    fmpz_set_si(t, y);
    status = gr_div_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int gr_generic_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    fmpq_t t;
    int status;

    if (fmpq_is_zero(y))   /* for non-rings? */
    {
        gr_ptr t;
        status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        status |= gr_set_fmpq(t, y, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(res, x, t, ctx);

        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    fmpq_init(t);
    fmpq_inv(t, y);
    status = gr_mul_fmpq(res, x, t, ctx);
    fmpq_clear(t);
    return status;
}

int gr_generic_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, y, y_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, x, t, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, x, x_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, t, y, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int gr_generic_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return gr_div(res, x, y, ctx);
}


/* Generic powering */

/* Assumes exp >= 2; res and tmp not not aliased with x. */
static int
gr_generic_pow_ui_binexp(gr_ptr res, gr_ptr tmp, gr_srcptr x, ulong exp, gr_ctx_t ctx)
{
    gr_ptr R, S, T;
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    int status;
    int zeros;
    ulong bit;

    status = GR_SUCCESS;

    /* Determine parity due to swaps */
    zeros = 0;
    bit = exp;
    while (bit > 1)
    {
        zeros += !(bit & 1);
        bit >>= 1;
    }

    if (zeros % 2)
    {
        R = res;
        S = tmp;
    }
    else
    {
        R = tmp;
        S = res;
    }

    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);

    status |= sqr(R, x, ctx);

    if (bit & exp)
    {
        status |= mul(S, R, x, ctx);
        T = R;
        R = S;
        S = T;
    }

    while (bit >>= 1)
    {
        status |= sqr(S, R, ctx);

        if (bit & exp)
        {
            status |= mul(R, S, x, ctx);
        }
        else
        {
            T = R;
            R = S;
            S = T;
        }
    }

    return status;
}

/* todo: optimize swaps */
static int
gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx)
{
    gr_ptr t, u;
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_swap_op swap = GR_SWAP_OP(ctx, SWAP);
    int status;
    slong i;

    status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);

    status |= gr_set(t, x, ctx);

    for (i = fmpz_bits(exp) - 2; i >= 0; i--)
    {
        status |= sqr(u, t, ctx);

        if (fmpz_tstbit(exp, i))
            status |= mul(t, u, x, ctx);
        else
            swap(t, u, ctx);
    }

    swap(res, t, ctx);

    GR_TMP_CLEAR2(t, u, ctx);

    return status;
}

int
gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
{
    int status;

    if (e == 0)
    {
        return gr_one(res, ctx);
    }
    else if (e == 1)
    {
        return gr_set(res, x, ctx);
    }
    else if (e == 2)
    {
        return gr_sqr(res, x, ctx);
    }
    else if (e == 4)
    {
        status = gr_sqr(res, x, ctx);
        status |= gr_sqr(res, res, ctx);
        return status;
    }
    else
    {
        gr_ptr t, u;

        if (res == x)
        {
            GR_TMP_INIT2(t, u, ctx);
            status = gr_set(u, x, ctx);
            status |= gr_generic_pow_ui_binexp(res, t, u, e, ctx);
            GR_TMP_CLEAR2(t, u, ctx);
        }
        else
        {
            GR_TMP_INIT(t, ctx);
            status = gr_generic_pow_ui_binexp(res, t, x, e, ctx);
            GR_TMP_CLEAR(t, ctx);
        }

        return status;
    }
}

/* todo: call gr_pow_ui instead of gr_generic_pow_ui? */
int
gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    if (e >= 0)
    {
        return gr_generic_pow_ui(res, x, e, ctx);
    }
    else
    {
        int status;

        /* todo: some heuristic for when we want to invert before/after powering */
        status = gr_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_generic_pow_ui(res, res, -e, ctx);

        return status;
    }
}

int
gr_generic_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t e, gr_ctx_t ctx)
{
    int status;

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);

        /* todo: some heuristic for when we want to invert before/after powering */
        status = gr_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_generic_pow_fmpz(res, res, f, ctx);

        fmpz_clear(f);
        return status;
    }

    if (*e == 0)
        return gr_one(res, ctx);
    else if (*e == 1)
        return gr_set(res, x, ctx);
    else if (*e == 2)
        return gr_sqr(res, x, ctx);
    else
        return gr_generic_pow_fmpz_binexp(res, x, e, ctx);
}

/* at least catch square roots -- todo: generalize to nth roots */
int
gr_generic_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(y)))
    {
        return gr_pow_fmpz(res, x, fmpq_numref(y), ctx);
    }
    else if (*fmpq_denref(y) == 2)
    {
        gr_ptr t;
        int status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        if (fmpz_sgn(fmpq_numref(y)) > 0)
        {
            status |= gr_sqrt(t, x, ctx);

            if (status == GR_SUCCESS)
                status |= gr_pow_fmpz(res, t, fmpq_numref(y), ctx);
        }
        else
        {
            status |= gr_rsqrt(t, x, ctx);

            if (status == GR_SUCCESS)
            {
                fmpz_t e;
                fmpz_init(e);
                fmpz_neg(e, fmpq_numref(y));
                status |= gr_pow_fmpz(res, t, e, ctx);
                fmpz_clear(e);
            }
        }

        if (status != GR_SUCCESS)
            status = GR_UNABLE;

        GR_TMP_CLEAR(t, ctx);
        return status;
    }
    else
    {
        gr_ptr t;
        int status;

        status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        status |= gr_set_fmpq(t, y, ctx);
        if (status == GR_SUCCESS)
            status = gr_pow(res, x, t, ctx);
        else
            status = GR_UNABLE;

        GR_TMP_CLEAR(t, ctx);
        return status;
    }
}

int gr_generic_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    if (y_ctx->which_ring == GR_CTX_FMPZ)
    {
        return gr_pow_fmpz(res, x, y, ctx);
    }
    else if (y_ctx->which_ring == GR_CTX_FMPQ)
    {
        return gr_pow_fmpq(res, x, y, ctx);
    }
    else
    {
        gr_ptr t;
        int status;

        status = GR_SUCCESS;

        GR_TMP_INIT(t, ctx);

        status |= gr_set_other(t, y, y_ctx, ctx);
        if (status == GR_SUCCESS)
            status = gr_pow(res, x, t, ctx);
        else
            status = GR_UNABLE;

        GR_TMP_CLEAR(t, ctx);
        return status;
    }
}

int gr_generic_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;

    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_set_other(t, x, x_ctx, ctx);
    if (status == GR_SUCCESS)
        status = gr_pow(res, t, y, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}


truth_t
gr_generic_is_square(gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return T_TRUE;

    if (gr_is_one(x, ctx) == T_TRUE)
        return T_TRUE;

    return GR_UNABLE;
}

int
gr_generic_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return GR_DOMAIN;

    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    {
        int status;

        status = gr_sqrt(res, x, ctx);
        if (status == GR_SUCCESS)
        {
            status = gr_inv(res, res, ctx);
            if (status == GR_SUCCESS)
                return status;
        }
    }

    return GR_UNABLE;
}

int
gr_generic_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return GR_UNABLE;
}

int
gr_generic_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    gr_ptr t, u;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, ctx);
    status |= gr_abs(t, x, ctx);
    status |= gr_abs(u, y, ctx);
    status |= gr_cmp(res, t, u, ctx);
    GR_TMP_CLEAR2(t, u, ctx);

    return status;
}

int
gr_generic_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_other(t, y, y_ctx, ctx);
    status |= gr_cmp(res, x, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_other(t, y, y_ctx, ctx);
    status |= gr_cmpabs(res, x, t, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = gr_exp(res, x, ctx);
    status |= gr_sub_ui(res, res, 1, ctx);
    return status;
}

int
gr_generic_exp2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_ui(t, 2, ctx);
    status |= gr_pow(res, t, x, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_exp10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_set_ui(t, 10, ctx);
    status |= gr_pow(res, t, x, ctx);
    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    int status = gr_add_ui(res, x, 1, ctx);
    status |= gr_log(res, res, ctx);
    return status;
}

int
gr_generic_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_one(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}

int
gr_generic_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
        return gr_zero(res, ctx);

    return GR_UNABLE;
}


int
gr_generic_bernoulli_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPQ)
    {
        if (0 && n <= 1000)
            bernoulli_cache_compute(n + 1);

        bernoulli_fmpq_ui(res, n);
        return GR_SUCCESS;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_bernoulli_ui(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else        /* compute via fmpq */
    {
        int status;
        fmpq_t t;
        fmpq_init(t);

        if (0 && n <= 1000)
            bernoulli_cache_compute(n + 1);
        bernoulli_fmpq_ui(t, n);

        status = gr_set_fmpq(res, t, ctx);
        fmpq_clear(t);
        return status;
    }
}

int
gr_generic_bernoulli_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*n) && *n >= 0)
    {
        return gr_bernoulli_ui(res, *n, ctx);
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_bernoulli_fmpz(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
gr_generic_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPQ)
    {
        if (0 && len <= 3000)
        {
            slong i;
            bernoulli_cache_compute(len);
            for (i = 0; i < len; i++)
                bernoulli_fmpq_ui(((fmpq *) res) + i, i);
        }
        else
        {
            /* todo: read as many as already cached ... */
            bernoulli_fmpq_vec_no_cache(res, 0, len);
        }

        return GR_SUCCESS;
    }

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute numerically via arb */
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        if (len > prec)
        {
            slong i, sz = ctx->sizeof_elem;
            int status = GR_SUCCESS;
            gr_ctx_t RR;
            arb_t t;

            gr_ctx_init_real_arb(RR, prec);

            arb_init(t);
            for (i = 0; i < len; i++)
            {
                arb_bernoulli_ui(t, i, prec);
                status |= gr_set_other(GR_ENTRY(res, i, sz), t, RR, ctx);
            }

            arb_clear(t);
            gr_ctx_clear(RR);

            return GR_SUCCESS;
        }
    }

    /* compute exactly via Q */
    {
        int status = GR_SUCCESS;
        slong i, sz = ctx->sizeof_elem;
        gr_ctx_t QQ;
        fmpq * t;
        gr_ctx_init_fmpq(QQ);
        GR_TMP_INIT_VEC(t, len, QQ);

        if (0 && len <= 3000)
        {
            slong i;
            bernoulli_cache_compute(len);
            for (i = 0; i < len; i++)
                bernoulli_fmpq_ui(t + i, i);
        }
        else
        {
            /* todo: read as many as already cached ... */
            bernoulli_fmpq_vec_no_cache(t, 0, len);
        }

        for (i = 0; i < len && status == GR_SUCCESS; i++)
            status |= gr_set_fmpq(GR_ENTRY(res, i, sz), t + i, ctx);

        GR_TMP_CLEAR_VEC(t, len, QQ);
        gr_ctx_clear(QQ);
        return status;
    }
}

void arb_fmpz_euler_number_ui(fmpz_t res, ulong n);

int
gr_generic_eulernum_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arb_fmpz_euler_number_ui(res, n);
        return GR_SUCCESS;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_eulernum_ui(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else        /* compute via fmpz */
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        arb_fmpz_euler_number_ui(t, n);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_eulernum_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*n) && *n >= 0)
    {
        return gr_eulernum_ui(res, *n, ctx);
    }
    else if (fmpz_sgn(n) < 0)
    {
        return GR_DOMAIN;
    }
    else if (fmpz_is_odd(n))
    {
        return gr_zero(res, ctx);
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute via arb */
    {
        gr_ctx_t RR;
        arb_t t;
        slong prec;
        int status = GR_SUCCESS;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));
        gr_ctx_init_real_arb(RR, prec);
        arb_init(t);
        status |= gr_eulernum_fmpz(t, n, RR);
        status |= gr_set_other(res, t, RR, ctx);
        arb_clear(t);
        gr_ctx_clear(RR);
        return status;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
gr_generic_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx)
{
    /* todo: arb based fmpz algorithm similar to bernoulli */
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_euler_number_vec(res, len);
        return GR_SUCCESS;
    }

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)  /* compute numerically via arb */
    {
        slong prec;
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&prec, ctx));

        if (len > prec)
        {
            slong i, sz = ctx->sizeof_elem;
            int status = GR_SUCCESS;
            gr_ctx_t RR;
            arb_t t;

            gr_ctx_init_real_arb(RR, prec);

            arb_init(t);
            for (i = 0; i < len; i++)
            {
                arb_euler_number_ui(t, i, prec);
                status |= gr_set_other(GR_ENTRY(res, i, sz), t, RR, ctx);
            }

            arb_clear(t);
            gr_ctx_clear(RR);

            return GR_SUCCESS;
        }
    }

    /* compute exactly via Z */
    {
        int status = GR_SUCCESS;
        slong i, sz = ctx->sizeof_elem;
        gr_ctx_t ZZ;
        fmpz * t;
        gr_ctx_init_fmpq(ZZ);
        GR_TMP_INIT_VEC(t, len, ZZ);

        arith_euler_number_vec(t, len);
        for (i = 0; i < len && status == GR_SUCCESS; i++)
            status |= gr_set_fmpz(GR_ENTRY(res, i, sz), t + i, ctx);

        GR_TMP_CLEAR_VEC(t, len, ZZ);
        gr_ctx_clear(ZZ);
        return status;
    }
}

int
gr_generic_stirling_s1u_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_stirling_number_1u(res, x, y);
        return GR_SUCCESS;
    }
    else
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        arith_stirling_number_1u(t, x, y);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_stirling_s1_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_stirling_number_1(res, x, y);
        return GR_SUCCESS;
    }
    else
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        arith_stirling_number_1(t, x, y);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_stirling_s2_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_stirling_number_2(res, x, y);
        return GR_SUCCESS;
    }
    else
    {
        int status;
        fmpz_t t;
        fmpz_init(t);
        arith_stirling_number_2(t, x, y);
        status = gr_set_fmpz(res, t, ctx);
        fmpz_clear(t);
        return status;
    }
}

int
gr_generic_stirling_s1u_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_stirling_number_1u_vec(res, x, len);
        return GR_SUCCESS;
    }
    else
    {
        int status = GR_SUCCESS;
        slong i, sz = ctx->sizeof_elem;
        gr_ctx_t ZZ;
        fmpz * t;
        gr_ctx_init_fmpz(ZZ);
        GR_TMP_INIT_VEC(t, len, ZZ);
        arith_stirling_number_1u_vec(t, x, len);
        /* todo: vector method */
        for (i = 0; i < len; i++)
            status |= gr_set_fmpz(GR_ENTRY(res, i, sz), t + i, ctx);
        GR_TMP_CLEAR_VEC(t, len, ZZ);
        gr_ctx_clear(ZZ);
        return status;
    }
}

int
gr_generic_stirling_s1_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_stirling_number_1_vec(res, x, len);
        return GR_SUCCESS;
    }
    else
    {
        int status = GR_SUCCESS;
        slong i, sz = ctx->sizeof_elem;
        gr_ctx_t ZZ;
        fmpz * t;
        gr_ctx_init_fmpz(ZZ);
        GR_TMP_INIT_VEC(t, len, ZZ);
        arith_stirling_number_1_vec(t, x, len);
        /* todo: vector method */
        for (i = 0; i < len; i++)
            status |= gr_set_fmpz(GR_ENTRY(res, i, sz), t + i, ctx);
        GR_TMP_CLEAR_VEC(t, len, ZZ);
        gr_ctx_clear(ZZ);
        return status;
    }
}

int
gr_generic_stirling_s2_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_FMPZ)
    {
        arith_stirling_number_2_vec(res, x, len);
        return GR_SUCCESS;
    }
    else
    {
        int status = GR_SUCCESS;
        slong i, sz = ctx->sizeof_elem;
        gr_ctx_t ZZ;
        fmpz * t;
        gr_ctx_init_fmpz(ZZ);
        GR_TMP_INIT_VEC(t, len, ZZ);
        arith_stirling_number_2_vec(t, x, len);
        /* todo: vector method */
        for (i = 0; i < len; i++)
            status |= gr_set_fmpz(GR_ENTRY(res, i, sz), t + i, ctx);
        GR_TMP_CLEAR_VEC(t, len, ZZ);
        gr_ctx_clear(ZZ);
        return status;
    }
}





/* Generic vector functions */

void
gr_generic_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT);
    slong i, sz;

    sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        init(GR_ENTRY(vec, i, sz), ctx);
}

void
gr_generic_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR);
    slong i, sz;

    sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        clear(GR_ENTRY(vec, i, sz), ctx);
}

void
gr_generic_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_swap_op swap = GR_SWAP_OP(ctx, SWAP);
    slong i, sz;

    sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        swap(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);
}

int
gr_generic_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_constant_op zero = GR_CONSTANT_OP(ctx, ZERO);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= zero(GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
gr_generic_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    gr_method_unary_op set = GR_UNARY_OP(ctx, SET);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= set(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

int
gr_generic_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    gr_method_unary_op neg = GR_UNARY_OP(ctx, NEG);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= neg(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

#define BINARY_OP(OP, op) \
int \
gr_generic_vec_ ## op(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) \
{ \
    gr_method_binary_op op = GR_BINARY_OP(ctx, OP); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), GR_ENTRY(src2, i, sz), ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_vec_ ## op ## _scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) \
{ \
    gr_method_binary_op op = GR_BINARY_OP(ctx, OP); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_vec_ ## op ## _scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) \
{ \
    gr_method_binary_op_si op = GR_BINARY_OP_SI(ctx, OP ## _SI); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_vec_ ## op ## _scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx) \
{ \
    gr_method_binary_op_ui op = GR_BINARY_OP_UI(ctx, OP ## _UI); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_vec_ ## op ## _scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx) \
{ \
    gr_method_binary_op_fmpz op = GR_BINARY_OP_FMPZ(ctx, OP ## _FMPZ); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_vec_ ## op ## _scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx) \
{ \
    gr_method_binary_op_fmpq op = GR_BINARY_OP_FMPQ(ctx, OP ## _FMPQ); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_scalar_ ## op ## _vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx) \
{ \
    gr_method_binary_op op = GR_BINARY_OP(ctx, OP); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(vec1, i, sz), c, GR_ENTRY(vec2, i, sz), ctx); \
 \
    return status; \
} \
int \
gr_generic_vec_ ## op ## _other(gr_ptr res, gr_srcptr src1, gr_srcptr src2, gr_ctx_t src2_ctx, slong len, gr_ctx_t ctx) \
{ \
    gr_method_binary_op_other op = GR_BINARY_OP_OTHER(ctx, OP ## _OTHER); \
    int status; \
    slong i, sz, sz2; \
 \
    sz = ctx->sizeof_elem; \
    sz2 = src2_ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), GR_ENTRY(src2, i, sz2), src2_ctx, ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_other_ ## op ## _vec(gr_ptr res, gr_srcptr src1, gr_ctx_t src1_ctx, gr_srcptr src2, slong len, gr_ctx_t ctx) \
{ \
    gr_method_other_binary_op op = GR_OTHER_BINARY_OP(ctx, OTHER_ ## OP); \
    int status; \
    slong i, sz, sz1; \
 \
    sz = ctx->sizeof_elem; \
    sz1 = src1_ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz1), src1_ctx, GR_ENTRY(src2, i, sz), ctx); \
 \
    return status; \
} \
 \
int \
gr_generic_vec_ ## op ## _scalar_other(gr_ptr res, gr_srcptr src1, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx) \
{ \
    gr_method_binary_op_other op = GR_BINARY_OP_OTHER(ctx, OP ## _OTHER); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), c, cctx, ctx); \
 \
    return status; \
} \
int \
gr_generic_scalar_other_ ## op ## _vec(gr_ptr res, gr_srcptr c, gr_ctx_t cctx, gr_srcptr src1, slong len, gr_ctx_t ctx) \
{ \
    gr_method_other_binary_op op = GR_OTHER_BINARY_OP(ctx, OTHER_ ## OP); \
    int status; \
    slong i, sz; \
 \
    sz = ctx->sizeof_elem; \
    status = GR_SUCCESS; \
 \
    for (i = 0; i < len; i++) \
        status |= op(GR_ENTRY(res, i, sz), c, cctx, GR_ENTRY(src1, i, sz), ctx); \
 \
    return status; \
} \


BINARY_OP(ADD, add)
BINARY_OP(SUB, sub)
BINARY_OP(MUL, mul)
BINARY_OP(DIV, div)
BINARY_OP(DIVEXACT, divexact)
BINARY_OP(POW, pow)

int
gr_generic_vec_mul_scalar_2exp_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    gr_method_binary_op_si op = GR_BINARY_OP_SI(ctx, MUL_2EXP_SI);
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= op(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), c, ctx);

    return status;
}

int
gr_generic_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= add(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_generic_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= sub(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_generic_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    gr_method_binary_op_si mul_si = GR_BINARY_OP_SI(ctx, MUL_SI);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul_si(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= add(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_generic_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    gr_method_binary_op_si mul_si = GR_BINARY_OP_SI(ctx, MUL_SI);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= mul_si(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= sub(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR(t, ctx);
    return status;
}

truth_t
gr_generic_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_predicate equal = GR_BINARY_PREDICATE(ctx, EQUAL);
    truth_t eq, this_eq;
    slong i, sz;

    sz = ctx->sizeof_elem;

    eq = T_TRUE;

    for (i = 0; i < len; i++)
    {
        this_eq = equal(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;

        if (this_eq == T_UNKNOWN)
            eq = T_UNKNOWN;
    }

    return eq;
}

int
gr_generic_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    truth_t eq, this_eq;
    slong i, sz;

    sz = ctx->sizeof_elem;

    eq = T_TRUE;

    for (i = 0; i < len; i++)
    {
        this_eq = is_zero(GR_ENTRY(vec, i, sz), ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;

        if (this_eq == T_UNKNOWN)
            eq = T_UNKNOWN;
    }

    return eq;
}

int
gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    if (initial == NULL)
    {
        status |= mul(res, vec1, vec2, ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul(t, vec1, vec2, ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_generic_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    if (initial == NULL)
    {
        status |= mul(res, vec1, GR_ENTRY(vec2, len - 1, sz), ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul(t, vec1, GR_ENTRY(vec2, len - 1, sz), ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, len - 1 - i, sz), ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_ui mul_ui = GR_BINARY_OP_UI(ctx, MUL_UI);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    if (initial == NULL)
    {
        status |= mul_ui(res, vec1, vec2[0], ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul_ui(t, vec1, vec2[0], ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul_ui(t, GR_ENTRY(vec1, i, sz), vec2[i], ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_generic_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_si mul_si = GR_BINARY_OP_SI(ctx, MUL_SI);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    if (initial == NULL)
    {
        status |= mul_si(res, vec1, vec2[0], ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul_si(t, vec1, vec2[0], ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul_si(t, GR_ENTRY(vec1, i, sz), vec2[i], ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

int
gr_generic_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op_fmpz mul_fmpz = GR_BINARY_OP_FMPZ(ctx, MUL_FMPZ);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    if (initial == NULL)
    {
        status |= mul_fmpz(res, vec1, vec2, ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= mul_fmpz(t, vec1, vec2, ctx);
        status |= add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= mul_fmpz(t, GR_ENTRY(vec1, i, sz), vec2 + i, ctx);
        status |= add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* todo: when to multiply repeatedly by x vs squaring? */
int
gr_generic_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            status |= gr_one(GR_ENTRY(res, i, sz), ctx);
        else if (i == 1)
            status |= gr_set(GR_ENTRY(res, i, sz), x, ctx);
        else if (i % 2 == 0)
            status |= sqr(GR_ENTRY(res, i, sz), GR_ENTRY(res, i / 2, sz), ctx);
        else
            status |= mul(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), x, ctx);
    }

    return status;
}

int
gr_generic_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;;

    for (i = 0; i < len && status == GR_SUCCESS; i++)
    {
        status |= gr_set_ui(GR_ENTRY(res, i, sz), i + 1, ctx);
        status |= gr_inv(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), ctx);
    }

    return status;
}


int
gr_generic_ctx_clear(gr_ctx_t ctx)
{
    return GR_SUCCESS;
}

/* Generic method implementations */
const gr_method_tab_input _gr_generic_methods[] =
{
    {GR_METHOD_CTX_CLEAR,               (gr_funcptr) gr_generic_ctx_clear},

    {GR_METHOD_CTX_IS_RING,             (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,      (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,        (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,         (gr_funcptr) gr_generic_ctx_predicate},

    {GR_METHOD_CTX_IS_ORDERED_RING,     (gr_funcptr) gr_generic_ctx_predicate},

    /* should be true for most rings */
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) gr_generic_ctx_predicate_true},

    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) gr_generic_ctx_predicate},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) gr_generic_ctx_predicate},

    {GR_METHOD_INIT,                    (gr_funcptr) gr_generic_init},
    {GR_METHOD_CLEAR,                   (gr_funcptr) gr_generic_clear},
    {GR_METHOD_SWAP,                    (gr_funcptr) gr_generic_swap},
    {GR_METHOD_SET_SHALLOW,             (gr_funcptr) gr_generic_set_shallow},
    {GR_METHOD_WRITE,                   (gr_funcptr) gr_generic_write},

    {GR_METHOD_RANDTEST,                (gr_funcptr) gr_generic_randtest},
    {GR_METHOD_RANDTEST_NOT_ZERO,       (gr_funcptr) gr_generic_randtest_not_zero},

    {GR_METHOD_ZERO,                    (gr_funcptr) gr_generic_zero},
    {GR_METHOD_ONE,                     (gr_funcptr) gr_generic_one},
    {GR_METHOD_NEG_ONE,                 (gr_funcptr) gr_generic_neg_one},

    {GR_METHOD_IS_ZERO,                 (gr_funcptr) gr_generic_is_zero},
    {GR_METHOD_IS_ONE,                  (gr_funcptr) gr_generic_is_one},
    {GR_METHOD_IS_NEG_ONE,              (gr_funcptr) gr_generic_is_neg_one},

    {GR_METHOD_EQUAL,                   (gr_funcptr) gr_generic_equal},

    {GR_METHOD_SET,                     (gr_funcptr) gr_generic_set},
    {GR_METHOD_SET_UI,                  (gr_funcptr) gr_generic_set_ui},
    {GR_METHOD_SET_SI,                  (gr_funcptr) gr_generic_set_si},
    {GR_METHOD_SET_FMPZ,                (gr_funcptr) gr_generic_set_fmpz},
    {GR_METHOD_SET_FMPQ,                (gr_funcptr) gr_generic_set_fmpq},

    {GR_METHOD_SET_OTHER,               (gr_funcptr) gr_generic_set_other},

    {GR_METHOD_NEG,                     (gr_funcptr) gr_generic_neg},

    {GR_METHOD_ADD,                     (gr_funcptr) gr_generic_add},
    {GR_METHOD_ADD_UI,                  (gr_funcptr) gr_generic_add_ui},
    {GR_METHOD_ADD_SI,                  (gr_funcptr) gr_generic_add_si},
    {GR_METHOD_ADD_FMPZ,                (gr_funcptr) gr_generic_add_fmpz},
    {GR_METHOD_ADD_FMPQ,                (gr_funcptr) gr_generic_add_fmpq},
    {GR_METHOD_ADD_OTHER,               (gr_funcptr) gr_generic_add_other},
    {GR_METHOD_OTHER_ADD,               (gr_funcptr) gr_generic_other_add},

    {GR_METHOD_SUB,                     (gr_funcptr) gr_generic_sub},
    {GR_METHOD_SUB_UI,                  (gr_funcptr) gr_generic_sub_ui},
    {GR_METHOD_SUB_SI,                  (gr_funcptr) gr_generic_sub_si},
    {GR_METHOD_SUB_FMPZ,                (gr_funcptr) gr_generic_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,                (gr_funcptr) gr_generic_sub_fmpq},
    {GR_METHOD_SUB_OTHER,               (gr_funcptr) gr_generic_sub_other},
    {GR_METHOD_OTHER_SUB,               (gr_funcptr) gr_generic_other_sub},

    {GR_METHOD_MUL,                     (gr_funcptr) gr_generic_mul},
    {GR_METHOD_MUL_UI,                  (gr_funcptr) gr_generic_mul_ui},
    {GR_METHOD_MUL_SI,                  (gr_funcptr) gr_generic_mul_si},
    {GR_METHOD_MUL_FMPZ,                (gr_funcptr) gr_generic_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,                (gr_funcptr) gr_generic_mul_fmpq},
    {GR_METHOD_MUL_OTHER,               (gr_funcptr) gr_generic_mul_other},
    {GR_METHOD_OTHER_MUL,               (gr_funcptr) gr_generic_other_mul},

    {GR_METHOD_ADDMUL,                  (gr_funcptr) gr_generic_addmul},
    {GR_METHOD_ADDMUL_UI,               (gr_funcptr) gr_generic_addmul_ui},
    {GR_METHOD_ADDMUL_SI,               (gr_funcptr) gr_generic_addmul_si},
    {GR_METHOD_ADDMUL_FMPZ,             (gr_funcptr) gr_generic_addmul_fmpz},
    {GR_METHOD_ADDMUL_FMPQ,             (gr_funcptr) gr_generic_addmul_fmpq},
    {GR_METHOD_ADDMUL_OTHER,            (gr_funcptr) gr_generic_addmul_other},

    {GR_METHOD_SUBMUL,                  (gr_funcptr) gr_generic_submul},
    {GR_METHOD_SUBMUL_UI,               (gr_funcptr) gr_generic_submul_ui},
    {GR_METHOD_SUBMUL_SI,               (gr_funcptr) gr_generic_submul_si},
    {GR_METHOD_SUBMUL_FMPZ,             (gr_funcptr) gr_generic_submul_fmpz},
    {GR_METHOD_SUBMUL_FMPQ,             (gr_funcptr) gr_generic_submul_fmpq},
    {GR_METHOD_SUBMUL_OTHER,            (gr_funcptr) gr_generic_submul_other},

    {GR_METHOD_MUL_TWO,                 (gr_funcptr) gr_generic_mul_two},
    {GR_METHOD_SQR,                     (gr_funcptr) gr_generic_sqr},

    {GR_METHOD_MUL_2EXP_SI,             (gr_funcptr) gr_generic_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,           (gr_funcptr) gr_generic_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ,      (gr_funcptr) gr_generic_set_fmpz_2exp_fmpz},
    {GR_METHOD_GET_FMPZ_2EXP_FMPZ,      (gr_funcptr) gr_generic_get_fmpz_2exp_fmpz},

    {GR_METHOD_DIV_UI,                  (gr_funcptr) gr_generic_div_ui},
    {GR_METHOD_DIV_SI,                  (gr_funcptr) gr_generic_div_si},
    {GR_METHOD_DIV_FMPZ,                (gr_funcptr) gr_generic_div_fmpz},
    {GR_METHOD_DIV_FMPQ,                (gr_funcptr) gr_generic_div_fmpq},
    {GR_METHOD_DIV_OTHER,               (gr_funcptr) gr_generic_div_other},
    {GR_METHOD_OTHER_DIV,               (gr_funcptr) gr_generic_other_div},

    {GR_METHOD_DIVEXACT,                (gr_funcptr) gr_generic_divexact},
    {GR_METHOD_DIVEXACT_UI,             (gr_funcptr) gr_generic_div_ui},
    {GR_METHOD_DIVEXACT_SI,             (gr_funcptr) gr_generic_div_si},
    {GR_METHOD_DIVEXACT_FMPZ,           (gr_funcptr) gr_generic_div_fmpz},
    {GR_METHOD_DIVEXACT_FMPQ,           (gr_funcptr) gr_generic_div_fmpq},
    {GR_METHOD_DIVEXACT_OTHER,          (gr_funcptr) gr_generic_div_other},
    {GR_METHOD_OTHER_DIVEXACT,          (gr_funcptr) gr_generic_other_div},

    {GR_METHOD_IS_INVERTIBLE,           (gr_funcptr) gr_generic_is_invertible},
    {GR_METHOD_INV,                     (gr_funcptr) gr_generic_inv},

    {GR_METHOD_POW_SI,                  (gr_funcptr) gr_generic_pow_si},
    {GR_METHOD_POW_UI,                  (gr_funcptr) gr_generic_pow_ui},
    {GR_METHOD_POW_FMPZ,                (gr_funcptr) gr_generic_pow_fmpz},
    {GR_METHOD_POW_OTHER,               (gr_funcptr) gr_generic_pow_other},
    {GR_METHOD_POW_FMPQ,                (gr_funcptr) gr_generic_pow_fmpq},
    {GR_METHOD_OTHER_POW,               (gr_funcptr) gr_generic_other_pow},

    {GR_METHOD_IS_SQUARE,               (gr_funcptr) gr_generic_is_square},
    {GR_METHOD_SQRT,                    (gr_funcptr) gr_generic_sqrt},
    {GR_METHOD_RSQRT,                   (gr_funcptr) gr_generic_rsqrt},

    {GR_METHOD_CMP,                     (gr_funcptr) gr_generic_cmp},
    {GR_METHOD_CMPABS,                  (gr_funcptr) gr_generic_cmpabs},
    {GR_METHOD_CMP_OTHER,               (gr_funcptr) gr_generic_cmp_other},
    {GR_METHOD_CMPABS_OTHER,            (gr_funcptr) gr_generic_cmpabs_other},

    {GR_METHOD_EXP,                     (gr_funcptr) gr_generic_exp},
    {GR_METHOD_EXPM1,                   (gr_funcptr) gr_generic_expm1},
    {GR_METHOD_EXP2,                    (gr_funcptr) gr_generic_exp2},
    {GR_METHOD_EXP10,                   (gr_funcptr) gr_generic_exp10},
    {GR_METHOD_LOG,                     (gr_funcptr) gr_generic_log},
    {GR_METHOD_LOG1P,                   (gr_funcptr) gr_generic_log1p},
    {GR_METHOD_SIN,                     (gr_funcptr) gr_generic_sin},
    {GR_METHOD_COS,                     (gr_funcptr) gr_generic_cos},
    {GR_METHOD_TAN,                     (gr_funcptr) gr_generic_tan},
    {GR_METHOD_ATAN,                    (gr_funcptr) gr_generic_atan},
    {GR_METHOD_ATANH,                   (gr_funcptr) gr_generic_atanh},

    {GR_METHOD_FAC,                     (gr_funcptr) gr_generic_fac},
    {GR_METHOD_FAC_UI,                  (gr_funcptr) gr_generic_fac_ui},
    {GR_METHOD_FAC_FMPZ,                (gr_funcptr) gr_generic_fac_fmpz},
    {GR_METHOD_FAC_VEC,                 (gr_funcptr) gr_generic_fac_vec},

    {GR_METHOD_RFAC,                    (gr_funcptr) gr_generic_rfac},
    {GR_METHOD_RFAC_UI,                 (gr_funcptr) gr_generic_rfac_ui},
    {GR_METHOD_RFAC_FMPZ,               (gr_funcptr) gr_generic_rfac_fmpz},
    {GR_METHOD_RFAC_VEC,                (gr_funcptr) gr_generic_rfac_vec},

    {GR_METHOD_RISING,                  (gr_funcptr) gr_generic_rising},
    {GR_METHOD_RISING_UI,               (gr_funcptr) gr_generic_rising_ui},
    {GR_METHOD_FALLING,                 (gr_funcptr) gr_generic_falling},
    {GR_METHOD_FALLING_UI,              (gr_funcptr) gr_generic_falling_ui},

    {GR_METHOD_BIN,                     (gr_funcptr) gr_generic_bin},
    {GR_METHOD_BIN_UI,                  (gr_funcptr) gr_generic_bin_ui},
    {GR_METHOD_BIN_UIUI,                (gr_funcptr) gr_generic_bin_uiui},
    {GR_METHOD_BIN_VEC,                 (gr_funcptr) gr_generic_bin_vec},
    {GR_METHOD_BIN_UI_VEC,              (gr_funcptr) gr_generic_bin_ui_vec},

    {GR_METHOD_DOUBLEFAC,               (gr_funcptr) gr_generic_doublefac},
    {GR_METHOD_DOUBLEFAC_UI,            (gr_funcptr) gr_generic_doublefac_ui},

    {GR_METHOD_HARMONIC,                (gr_funcptr) gr_generic_harmonic},
    {GR_METHOD_HARMONIC_UI,             (gr_funcptr) gr_generic_harmonic_ui},

    {GR_METHOD_BETA,                    (gr_funcptr) gr_generic_beta},

    {GR_METHOD_BERNOULLI_UI,            (gr_funcptr) gr_generic_bernoulli_ui},
    {GR_METHOD_BERNOULLI_FMPZ,          (gr_funcptr) gr_generic_bernoulli_fmpz},
    {GR_METHOD_BERNOULLI_VEC,           (gr_funcptr) gr_generic_bernoulli_vec},

    {GR_METHOD_EULERNUM_UI,             (gr_funcptr) gr_generic_eulernum_ui},
    {GR_METHOD_EULERNUM_FMPZ,           (gr_funcptr) gr_generic_eulernum_fmpz},
    {GR_METHOD_EULERNUM_VEC,            (gr_funcptr) gr_generic_eulernum_vec},

    {GR_METHOD_FIB_UI,                  (gr_funcptr) gr_generic_fib_ui},
    {GR_METHOD_FIB_FMPZ,                (gr_funcptr) gr_generic_fib_fmpz},
    {GR_METHOD_FIB_VEC,                 (gr_funcptr) gr_generic_fib_vec},

    {GR_METHOD_BELLNUM_UI,              (gr_funcptr) gr_generic_bellnum_ui},
    {GR_METHOD_BELLNUM_FMPZ,            (gr_funcptr) gr_generic_bellnum_fmpz},
    {GR_METHOD_BELLNUM_VEC,             (gr_funcptr) gr_generic_bellnum_vec},

    {GR_METHOD_PARTITIONS_UI,           (gr_funcptr) gr_generic_partitions_ui},
    {GR_METHOD_PARTITIONS_FMPZ,         (gr_funcptr) gr_generic_partitions_fmpz},
    {GR_METHOD_PARTITIONS_VEC,          (gr_funcptr) gr_generic_partitions_vec},

    {GR_METHOD_STIRLING_S1U_UIUI,       (gr_funcptr) gr_generic_stirling_s1u_uiui},
    {GR_METHOD_STIRLING_S1_UIUI,        (gr_funcptr) gr_generic_stirling_s1_uiui},
    {GR_METHOD_STIRLING_S2_UIUI,        (gr_funcptr) gr_generic_stirling_s2_uiui},
    {GR_METHOD_STIRLING_S1U_UI_VEC,     (gr_funcptr) gr_generic_stirling_s1u_ui_vec},
    {GR_METHOD_STIRLING_S1_UI_VEC,      (gr_funcptr) gr_generic_stirling_s1_ui_vec},
    {GR_METHOD_STIRLING_S2_UI_VEC,      (gr_funcptr) gr_generic_stirling_s2_ui_vec},

    {GR_METHOD_CHEBYSHEV_T_FMPZ,        (gr_funcptr) gr_generic_chebyshev_t_fmpz},
    {GR_METHOD_CHEBYSHEV_U_FMPZ,        (gr_funcptr) gr_generic_chebyshev_u_fmpz},
 
    {GR_METHOD_HILBERT_CLASS_POLY,      (gr_funcptr) gr_generic_hilbert_class_poly},

    {GR_METHOD_VEC_INIT,                (gr_funcptr) gr_generic_vec_init},
    {GR_METHOD_VEC_CLEAR,               (gr_funcptr) gr_generic_vec_clear},
    {GR_METHOD_VEC_SWAP,                (gr_funcptr) gr_generic_vec_swap},
    {GR_METHOD_VEC_ZERO,                (gr_funcptr) gr_generic_vec_zero},
    {GR_METHOD_VEC_SET,                 (gr_funcptr) gr_generic_vec_set},
    {GR_METHOD_VEC_NEG,                 (gr_funcptr) gr_generic_vec_neg},

    {GR_METHOD_VEC_ADD,                 (gr_funcptr) gr_generic_vec_add},
    {GR_METHOD_VEC_SUB,                 (gr_funcptr) gr_generic_vec_sub},
    {GR_METHOD_VEC_MUL,                 (gr_funcptr) gr_generic_vec_mul},
    {GR_METHOD_VEC_DIV,                 (gr_funcptr) gr_generic_vec_div},
    {GR_METHOD_VEC_DIVEXACT,            (gr_funcptr) gr_generic_vec_divexact},
    {GR_METHOD_VEC_POW,                 (gr_funcptr) gr_generic_vec_pow},

    {GR_METHOD_VEC_ADD_SCALAR,          (gr_funcptr) gr_generic_vec_add_scalar},
    {GR_METHOD_VEC_SUB_SCALAR,          (gr_funcptr) gr_generic_vec_sub_scalar},
    {GR_METHOD_VEC_MUL_SCALAR,          (gr_funcptr) gr_generic_vec_mul_scalar},
    {GR_METHOD_VEC_DIV_SCALAR,          (gr_funcptr) gr_generic_vec_div_scalar},
    {GR_METHOD_VEC_DIVEXACT_SCALAR,     (gr_funcptr) gr_generic_vec_divexact_scalar},
    {GR_METHOD_VEC_POW_SCALAR,          (gr_funcptr) gr_generic_vec_pow_scalar},

    {GR_METHOD_VEC_ADD_SCALAR_SI,       (gr_funcptr) gr_generic_vec_add_scalar_si},
    {GR_METHOD_VEC_SUB_SCALAR_SI,       (gr_funcptr) gr_generic_vec_sub_scalar_si},
    {GR_METHOD_VEC_MUL_SCALAR_SI,       (gr_funcptr) gr_generic_vec_mul_scalar_si},
    {GR_METHOD_VEC_DIV_SCALAR_SI,       (gr_funcptr) gr_generic_vec_div_scalar_si},
    {GR_METHOD_VEC_DIVEXACT_SCALAR_SI,  (gr_funcptr) gr_generic_vec_divexact_scalar_si},
    {GR_METHOD_VEC_POW_SCALAR_SI,       (gr_funcptr) gr_generic_vec_pow_scalar_si},

    {GR_METHOD_VEC_ADD_SCALAR_UI,       (gr_funcptr) gr_generic_vec_add_scalar_ui},
    {GR_METHOD_VEC_SUB_SCALAR_UI,       (gr_funcptr) gr_generic_vec_sub_scalar_ui},
    {GR_METHOD_VEC_MUL_SCALAR_UI,       (gr_funcptr) gr_generic_vec_mul_scalar_ui},
    {GR_METHOD_VEC_DIV_SCALAR_UI,       (gr_funcptr) gr_generic_vec_div_scalar_ui},
    {GR_METHOD_VEC_DIVEXACT_SCALAR_UI,  (gr_funcptr) gr_generic_vec_divexact_scalar_ui},
    {GR_METHOD_VEC_POW_SCALAR_UI,       (gr_funcptr) gr_generic_vec_pow_scalar_ui},

    {GR_METHOD_VEC_ADD_SCALAR_FMPZ,     (gr_funcptr) gr_generic_vec_add_scalar_fmpz},
    {GR_METHOD_VEC_SUB_SCALAR_FMPZ,     (gr_funcptr) gr_generic_vec_sub_scalar_fmpz},
    {GR_METHOD_VEC_MUL_SCALAR_FMPZ,     (gr_funcptr) gr_generic_vec_mul_scalar_fmpz},
    {GR_METHOD_VEC_DIV_SCALAR_FMPZ,     (gr_funcptr) gr_generic_vec_div_scalar_fmpz},
    {GR_METHOD_VEC_DIVEXACT_SCALAR_FMPZ,(gr_funcptr) gr_generic_vec_divexact_scalar_fmpz},
    {GR_METHOD_VEC_POW_SCALAR_FMPZ,     (gr_funcptr) gr_generic_vec_pow_scalar_fmpz},

    {GR_METHOD_VEC_ADD_SCALAR_FMPQ,     (gr_funcptr) gr_generic_vec_add_scalar_fmpq},
    {GR_METHOD_VEC_SUB_SCALAR_FMPQ,     (gr_funcptr) gr_generic_vec_sub_scalar_fmpq},
    {GR_METHOD_VEC_MUL_SCALAR_FMPQ,     (gr_funcptr) gr_generic_vec_mul_scalar_fmpq},
    {GR_METHOD_VEC_DIV_SCALAR_FMPQ,     (gr_funcptr) gr_generic_vec_div_scalar_fmpq},
    {GR_METHOD_VEC_DIVEXACT_SCALAR_FMPQ,(gr_funcptr) gr_generic_vec_divexact_scalar_fmpq},
    {GR_METHOD_VEC_POW_SCALAR_FMPQ,     (gr_funcptr) gr_generic_vec_pow_scalar_fmpq},

    {GR_METHOD_SCALAR_ADD_VEC,          (gr_funcptr) gr_generic_scalar_add_vec},
    {GR_METHOD_SCALAR_SUB_VEC,          (gr_funcptr) gr_generic_scalar_sub_vec},
    {GR_METHOD_SCALAR_MUL_VEC,          (gr_funcptr) gr_generic_scalar_mul_vec},
    {GR_METHOD_SCALAR_DIV_VEC,          (gr_funcptr) gr_generic_scalar_div_vec},
    {GR_METHOD_SCALAR_DIVEXACT_VEC,     (gr_funcptr) gr_generic_scalar_divexact_vec},
    {GR_METHOD_SCALAR_POW_VEC,          (gr_funcptr) gr_generic_scalar_pow_vec},

    {GR_METHOD_VEC_ADD_OTHER,          (gr_funcptr) gr_generic_vec_add_other},
    {GR_METHOD_VEC_SUB_OTHER,          (gr_funcptr) gr_generic_vec_sub_other},
    {GR_METHOD_VEC_MUL_OTHER,          (gr_funcptr) gr_generic_vec_mul_other},
    {GR_METHOD_VEC_DIV_OTHER,          (gr_funcptr) gr_generic_vec_div_other},
    {GR_METHOD_VEC_DIVEXACT_OTHER,     (gr_funcptr) gr_generic_vec_divexact_other},
    {GR_METHOD_VEC_POW_OTHER,          (gr_funcptr) gr_generic_vec_pow_other},

    {GR_METHOD_OTHER_ADD_VEC,          (gr_funcptr) gr_generic_other_add_vec},
    {GR_METHOD_OTHER_SUB_VEC,          (gr_funcptr) gr_generic_other_sub_vec},
    {GR_METHOD_OTHER_MUL_VEC,          (gr_funcptr) gr_generic_other_mul_vec},
    {GR_METHOD_OTHER_DIV_VEC,          (gr_funcptr) gr_generic_other_div_vec},
    {GR_METHOD_OTHER_DIVEXACT_VEC,     (gr_funcptr) gr_generic_other_divexact_vec},
    {GR_METHOD_OTHER_POW_VEC,          (gr_funcptr) gr_generic_other_pow_vec},

    {GR_METHOD_VEC_ADD_SCALAR_OTHER,          (gr_funcptr) gr_generic_vec_add_scalar_other},
    {GR_METHOD_VEC_SUB_SCALAR_OTHER,          (gr_funcptr) gr_generic_vec_sub_scalar_other},
    {GR_METHOD_VEC_MUL_SCALAR_OTHER,          (gr_funcptr) gr_generic_vec_mul_scalar_other},
    {GR_METHOD_VEC_DIV_SCALAR_OTHER,          (gr_funcptr) gr_generic_vec_div_scalar_other},
    {GR_METHOD_VEC_DIVEXACT_SCALAR_OTHER,     (gr_funcptr) gr_generic_vec_divexact_scalar_other},
    {GR_METHOD_VEC_POW_SCALAR_OTHER,          (gr_funcptr) gr_generic_vec_pow_scalar_other},

    {GR_METHOD_SCALAR_OTHER_ADD_VEC,          (gr_funcptr) gr_generic_scalar_other_add_vec},
    {GR_METHOD_SCALAR_OTHER_SUB_VEC,          (gr_funcptr) gr_generic_scalar_other_sub_vec},
    {GR_METHOD_SCALAR_OTHER_MUL_VEC,          (gr_funcptr) gr_generic_scalar_other_mul_vec},
    {GR_METHOD_SCALAR_OTHER_DIV_VEC,          (gr_funcptr) gr_generic_scalar_other_div_vec},
    {GR_METHOD_SCALAR_OTHER_DIVEXACT_VEC,     (gr_funcptr) gr_generic_scalar_other_divexact_vec},

    {GR_METHOD_VEC_MUL_SCALAR_2EXP_SI,       (gr_funcptr) gr_generic_vec_mul_scalar_2exp_si},

    {GR_METHOD_VEC_ADDMUL_SCALAR,       (gr_funcptr) gr_generic_vec_scalar_addmul},
    {GR_METHOD_VEC_SUBMUL_SCALAR,       (gr_funcptr) gr_generic_vec_scalar_submul},
    {GR_METHOD_VEC_ADDMUL_SCALAR_SI,    (gr_funcptr) gr_generic_vec_scalar_addmul_si},
    {GR_METHOD_VEC_SUBMUL_SCALAR_SI,    (gr_funcptr) gr_generic_vec_scalar_submul_si},

    {GR_METHOD_VEC_EQUAL,               (gr_funcptr) gr_generic_vec_equal},
    {GR_METHOD_VEC_IS_ZERO,             (gr_funcptr) gr_generic_vec_is_zero},

    {GR_METHOD_VEC_SUM,                 (gr_funcptr) _gr_vec_sum_generic},
    {GR_METHOD_VEC_PRODUCT,             (gr_funcptr) _gr_vec_product_generic},

    {GR_METHOD_VEC_DOT,                 (gr_funcptr) gr_generic_vec_dot},
    {GR_METHOD_VEC_DOT_REV,             (gr_funcptr) gr_generic_vec_dot_rev},
    {GR_METHOD_VEC_DOT_UI,              (gr_funcptr) gr_generic_vec_dot_ui},
    {GR_METHOD_VEC_DOT_SI,              (gr_funcptr) gr_generic_vec_dot_si},
    {GR_METHOD_VEC_DOT_FMPZ,            (gr_funcptr) gr_generic_vec_dot_fmpz},

    {GR_METHOD_VEC_RECIPROCALS,         (gr_funcptr) gr_generic_vec_reciprocals},
    {GR_METHOD_VEC_SET_POWERS,          (gr_funcptr) gr_generic_vec_set_powers},

    {GR_METHOD_POLY_MULLOW,             (gr_funcptr) _gr_poly_mullow_generic},
    {GR_METHOD_POLY_DIVREM,             (gr_funcptr) _gr_poly_divrem_generic},
    {GR_METHOD_POLY_INV_SERIES,         (gr_funcptr) _gr_poly_inv_series_generic},
    {GR_METHOD_POLY_DIV_SERIES,         (gr_funcptr) _gr_poly_div_series_generic},
    {GR_METHOD_POLY_RSQRT_SERIES,       (gr_funcptr) _gr_poly_rsqrt_series_generic},
    {GR_METHOD_POLY_SQRT_SERIES,        (gr_funcptr) _gr_poly_sqrt_series_generic},
    {GR_METHOD_POLY_EXP_SERIES,         (gr_funcptr) _gr_poly_exp_series_generic},

    {GR_METHOD_MAT_MUL,                 (gr_funcptr) gr_mat_mul_generic},
    {GR_METHOD_MAT_DET,                 (gr_funcptr) gr_mat_det_generic},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,  (gr_funcptr) gr_mat_find_nonzero_pivot_generic},
    {GR_METHOD_MAT_DIAGONALIZATION,     (gr_funcptr) gr_mat_diagonalization_generic},

    {0,                                 (gr_funcptr) NULL},
};

void
gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab)
{
    slong i;

    for (i = 0; i < GR_METHOD_TAB_SIZE; i++)
        methods[i] = (gr_funcptr) gr_not_implemented;

    /* Assign generic methods as fallbacks */
    for (i = 0; ; i++)
    {
        if (_gr_generic_methods[i].function == NULL)
            break;

        if (_gr_generic_methods[i].index >= GR_METHOD_TAB_SIZE)
            abort();

        methods[_gr_generic_methods[i].index] = _gr_generic_methods[i].function;
    }

    for (i = 0; ; i++)
    {
        if (tab[i].function == NULL)
            break;

        if (tab[i].index >= GR_METHOD_TAB_SIZE)
            abort();

        methods[tab[i].index] = tab[i].function;
    }
}
