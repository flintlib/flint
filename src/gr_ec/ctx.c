/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr.h"
#include "gr_ec.h"

static void
_gr_ec_ctx_alloc(gr_ec_ctx_t ctx, gr_ctx_t base_ring)
{
    ctx->base_ring = base_ring;
    ctx->coeffs = flint_malloc(GR_EC_CTX_NUM_COEFFS * base_ring->sizeof_elem);
    _gr_vec_init(ctx->coeffs, GR_EC_CTX_NUM_COEFFS, base_ring);
    ctx->model = GR_EC_LONG_WEIERSTRASS;
}

void
gr_ec_ctx_clear(gr_ec_ctx_t ctx)
{
    _gr_vec_clear(ctx->coeffs, GR_EC_CTX_NUM_COEFFS, GR_EC_ELEM_CTX(ctx));
    flint_free(ctx->coeffs);
    ctx->coeffs = NULL;
}

/* b-invariants, discriminant and model, from the a-invariants */
static int
_gr_ec_ctx_derive(gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_srcptr a1 = GR_EC_A1(ctx), a2 = GR_EC_A2(ctx), a3 = GR_EC_A3(ctx);
    gr_srcptr a4 = GR_EC_A4(ctx), a6 = GR_EC_A6(ctx);
    gr_ptr b2 = GR_EC_B2(ctx), b4 = GR_EC_B4(ctx), b6 = GR_EC_B6(ctx);
    gr_ptr b8 = GR_EC_B8(ctx), disc = GR_EC_DISC(ctx);
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, R);

    /* b2 = a1^2 + 4 a2 */
    status |= gr_sqr(b2, a1, R);
    status |= gr_mul_ui(t, a2, 4, R);
    status |= gr_add(b2, b2, t, R);

    /* b4 = 2 a4 + a1 a3 */
    status |= gr_mul_two(b4, a4, R);
    status |= gr_mul(t, a1, a3, R);
    status |= gr_add(b4, b4, t, R);

    /* b6 = a3^2 + 4 a6 */
    status |= gr_sqr(b6, a3, R);
    status |= gr_mul_ui(t, a6, 4, R);
    status |= gr_add(b6, b6, t, R);

    /* b8 = a1^2 a6 + 4 a2 a6 - a1 a3 a4 + a2 a3^2 - a4^2 */
    status |= gr_sqr(t, a1, R);
    status |= gr_mul(b8, t, a6, R);
    status |= gr_mul(t, a2, a6, R);
    status |= gr_mul_ui(t, t, 4, R);
    status |= gr_add(b8, b8, t, R);
    status |= gr_mul(t, a1, a3, R);
    status |= gr_mul(t, t, a4, R);
    status |= gr_sub(b8, b8, t, R);
    status |= gr_sqr(t, a3, R);
    status |= gr_mul(t, t, a2, R);
    status |= gr_add(b8, b8, t, R);
    status |= gr_sqr(t, a4, R);
    status |= gr_sub(b8, b8, t, R);

    /* disc = -b2^2 b8 - 8 b4^3 - 27 b6^2 + 9 b2 b4 b6 */
    status |= gr_sqr(t, b2, R);
    status |= gr_mul(t, t, b8, R);
    status |= gr_neg(disc, t, R);
    status |= gr_sqr(t, b4, R);
    status |= gr_mul(t, t, b4, R);
    status |= gr_mul_ui(t, t, 8, R);
    status |= gr_sub(disc, disc, t, R);
    status |= gr_sqr(t, b6, R);
    status |= gr_mul_ui(t, t, 27, R);
    status |= gr_sub(disc, disc, t, R);
    status |= gr_mul(t, b2, b4, R);
    status |= gr_mul(t, t, b6, R);
    status |= gr_mul_ui(t, t, 9, R);
    status |= gr_add(disc, disc, t, R);

    GR_TMP_CLEAR(t, R);

    if (status == GR_SUCCESS)
    {
        if (gr_is_zero(a1, R) == T_TRUE && gr_is_zero(a2, R) == T_TRUE
                && gr_is_zero(a3, R) == T_TRUE)
            ctx->model = GR_EC_SHORT_WEIERSTRASS;
        else
            ctx->model = GR_EC_LONG_WEIERSTRASS;
    }

    return status;
}

/* finish initialization; clears ctx and returns the flag on failure */
static int
_gr_ec_ctx_finish(gr_ec_ctx_t ctx, int status)
{
    if (status == GR_SUCCESS)
        status = _gr_ec_ctx_derive(ctx);

    if (status == GR_SUCCESS
            && gr_is_zero(GR_EC_DISC(ctx), GR_EC_ELEM_CTX(ctx)) == T_TRUE)
        status = GR_DOMAIN;

    if (status != GR_SUCCESS)
        gr_ec_ctx_clear(ctx);

    return status;
}

int
gr_ec_ctx_init(gr_ec_ctx_t ctx, gr_ctx_t base_ring, gr_srcptr a1, gr_srcptr a2,
        gr_srcptr a3, gr_srcptr a4, gr_srcptr a6)
{
    int status = GR_SUCCESS;

    if (gr_ctx_is_commutative_ring(base_ring) == T_FALSE)
        return GR_DOMAIN;

    _gr_ec_ctx_alloc(ctx, base_ring);

    status |= gr_set(GR_EC_A1(ctx), a1, base_ring);
    status |= gr_set(GR_EC_A2(ctx), a2, base_ring);
    status |= gr_set(GR_EC_A3(ctx), a3, base_ring);
    status |= gr_set(GR_EC_A4(ctx), a4, base_ring);
    status |= gr_set(GR_EC_A6(ctx), a6, base_ring);

    return _gr_ec_ctx_finish(ctx, status);
}

int
gr_ec_ctx_init_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring, slong a1, slong a2,
        slong a3, slong a4, slong a6)
{
    int status = GR_SUCCESS;

    if (gr_ctx_is_commutative_ring(base_ring) == T_FALSE)
        return GR_DOMAIN;

    _gr_ec_ctx_alloc(ctx, base_ring);

    status |= gr_set_si(GR_EC_A1(ctx), a1, base_ring);
    status |= gr_set_si(GR_EC_A2(ctx), a2, base_ring);
    status |= gr_set_si(GR_EC_A3(ctx), a3, base_ring);
    status |= gr_set_si(GR_EC_A4(ctx), a4, base_ring);
    status |= gr_set_si(GR_EC_A6(ctx), a6, base_ring);

    return _gr_ec_ctx_finish(ctx, status);
}

int
gr_ec_ctx_init_short_weierstrass(gr_ec_ctx_t ctx, gr_ctx_t base_ring,
        gr_srcptr a4, gr_srcptr a6)
{
    int status = GR_SUCCESS;

    if (gr_ctx_is_commutative_ring(base_ring) == T_FALSE)
        return GR_DOMAIN;

    _gr_ec_ctx_alloc(ctx, base_ring);

    status |= gr_zero(GR_EC_A1(ctx), base_ring);
    status |= gr_zero(GR_EC_A2(ctx), base_ring);
    status |= gr_zero(GR_EC_A3(ctx), base_ring);
    status |= gr_set(GR_EC_A4(ctx), a4, base_ring);
    status |= gr_set(GR_EC_A6(ctx), a6, base_ring);

    return _gr_ec_ctx_finish(ctx, status);
}

int
gr_ec_ctx_init_short_weierstrass_si(gr_ec_ctx_t ctx, gr_ctx_t base_ring,
        slong a4, slong a6)
{
    int status = GR_SUCCESS;

    if (gr_ctx_is_commutative_ring(base_ring) == T_FALSE)
        return GR_DOMAIN;

    _gr_ec_ctx_alloc(ctx, base_ring);

    status |= gr_zero(GR_EC_A1(ctx), base_ring);
    status |= gr_zero(GR_EC_A2(ctx), base_ring);
    status |= gr_zero(GR_EC_A3(ctx), base_ring);
    status |= gr_set_si(GR_EC_A4(ctx), a4, base_ring);
    status |= gr_set_si(GR_EC_A6(ctx), a6, base_ring);

    return _gr_ec_ctx_finish(ctx, status);
}

int
gr_ec_ctx_init_randtest(gr_ec_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring)
{
    slong iter;

    if (gr_ctx_is_commutative_ring(base_ring) == T_FALSE)
        return GR_DOMAIN;

    for (iter = 0; iter < 20; iter++)
    {
        int status = GR_SUCCESS;
        int shrt = n_randint(state, 2);

        _gr_ec_ctx_alloc(ctx, base_ring);

        if (shrt)
        {
            status |= gr_zero(GR_EC_A1(ctx), base_ring);
            status |= gr_zero(GR_EC_A2(ctx), base_ring);
            status |= gr_zero(GR_EC_A3(ctx), base_ring);
        }
        else
        {
            status |= gr_randtest(GR_EC_A1(ctx), state, base_ring);
            status |= gr_randtest(GR_EC_A2(ctx), state, base_ring);
            status |= gr_randtest(GR_EC_A3(ctx), state, base_ring);
        }

        status |= gr_randtest(GR_EC_A4(ctx), state, base_ring);
        status |= gr_randtest(GR_EC_A6(ctx), state, base_ring);

        status = _gr_ec_ctx_finish(ctx, status);

        if (status == GR_SUCCESS)
            return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
gr_ec_ctx_write(gr_stream_t out, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;
    slong i;

    status |= gr_stream_write(out, "Elliptic curve with a-invariants [");

    for (i = 0; i < 5; i++)
    {
        if (i != 0)
            status |= gr_stream_write(out, ", ");
        status |= gr_write(out, GR_EC_COEFF(ctx, i), R);
    }

    status |= gr_stream_write(out, "] over ");
    status |= gr_ctx_write(out, R);

    return status;
}

int
gr_ec_ctx_get_str(char ** res, gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    int status;
    gr_stream_init_str(out);
    status = gr_ec_ctx_write(out, ctx);
    *res = out->s;
    return status;
}

int
gr_ec_ctx_print(gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_ec_ctx_write(out, ctx);
}

int
gr_ec_ctx_a_invariants(gr_ptr a1, gr_ptr a2, gr_ptr a3, gr_ptr a4, gr_ptr a6,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(a1, GR_EC_A1(ctx), R);
    status |= gr_set(a2, GR_EC_A2(ctx), R);
    status |= gr_set(a3, GR_EC_A3(ctx), R);
    status |= gr_set(a4, GR_EC_A4(ctx), R);
    status |= gr_set(a6, GR_EC_A6(ctx), R);

    return status;
}

int
gr_ec_ctx_b_invariants(gr_ptr b2, gr_ptr b4, gr_ptr b6, gr_ptr b8,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(b2, GR_EC_B2(ctx), R);
    status |= gr_set(b4, GR_EC_B4(ctx), R);
    status |= gr_set(b6, GR_EC_B6(ctx), R);
    status |= gr_set(b8, GR_EC_B8(ctx), R);

    return status;
}

int
gr_ec_ctx_c_invariants(gr_ptr c4, gr_ptr c6, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_srcptr b2 = GR_EC_B2(ctx), b4 = GR_EC_B4(ctx), b6 = GR_EC_B6(ctx);
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, R);

    /* c4 = b2^2 - 24 b4 */
    status |= gr_sqr(c4, b2, R);
    status |= gr_mul_ui(t, b4, 24, R);
    status |= gr_sub(c4, c4, t, R);

    /* c6 = -b2^3 + 36 b2 b4 - 216 b6 */
    status |= gr_sqr(t, b2, R);
    status |= gr_mul(t, t, b2, R);
    status |= gr_neg(c6, t, R);
    status |= gr_mul(t, b2, b4, R);
    status |= gr_mul_ui(t, t, 36, R);
    status |= gr_add(c6, c6, t, R);
    status |= gr_mul_ui(t, b6, 216, R);
    status |= gr_sub(c6, c6, t, R);

    GR_TMP_CLEAR(t, R);

    return status;
}

int
gr_ec_ctx_discriminant(gr_ptr res, gr_ec_ctx_t ctx)
{
    return gr_set(res, GR_EC_DISC(ctx), GR_EC_ELEM_CTX(ctx));
}

int
gr_ec_ctx_j_invariant(gr_ptr res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr c4, c6;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(c4, c6, R);

    status |= gr_ec_ctx_c_invariants(c4, c6, ctx);

    /* j = c4^3 / disc */
    status |= gr_sqr(c6, c4, R);
    status |= gr_mul(c6, c6, c4, R);

    if (status == GR_SUCCESS)
    {
        truth_t inv = gr_is_invertible(GR_EC_DISC(ctx), R);

        if (inv == T_FALSE)
            status = GR_DOMAIN;
        else
            status = gr_div(res, c6, GR_EC_DISC(ctx), R);
    }

    GR_TMP_CLEAR2(c4, c6, R);

    return status;
}

truth_t
gr_ec_ctx_is_smooth(gr_ec_ctx_t ctx)
{
    return truth_not(gr_is_zero(GR_EC_DISC(ctx), GR_EC_ELEM_CTX(ctx)));
}
