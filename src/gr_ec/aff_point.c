/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "ulong_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_ec.h"

#define AX(P) GR_EC_AFF_POINT_X(P, ctx)
#define AY(P) GR_EC_AFF_POINT_Y(P, ctx)

#define T(i) GR_ENTRY(t, i, sz)

void
gr_ec_aff_point_init(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);

    P->coords = flint_malloc(2 * R->sizeof_elem);
    _gr_vec_init(P->coords, 2, R);
    P->is_infinity = T_TRUE;
}

void
gr_ec_aff_point_clear(gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
{
    _gr_vec_clear(P->coords, 2, GR_EC_ELEM_CTX(ctx));
    flint_free(P->coords);
    P->coords = NULL;
}

int
gr_ec_aff_point_set(gr_ec_aff_point_t res, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (P->is_infinity != T_TRUE)
        status = _gr_vec_set(res->coords, P->coords, 2, GR_EC_ELEM_CTX(ctx));

    res->is_infinity = P->is_infinity;

    return (P->is_infinity == T_UNKNOWN) ? GR_UNABLE : status;
}

int
gr_ec_aff_point_zero(gr_ec_aff_point_t res, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    res->is_infinity = T_TRUE;
    return GR_SUCCESS;
}

int
_gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(AX(res), x, R);
    status |= gr_set(AY(res), y, R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

int
gr_ec_aff_point_set_affine(gr_ec_aff_point_t res, gr_srcptr x, gr_srcptr y,
        gr_ec_ctx_t ctx)
{
    int status = _gr_ec_aff_point_set_affine(res, x, y, ctx);

    if (status == GR_SUCCESS)
    {
        truth_t on = gr_ec_aff_point_is_on_curve(res, ctx);

        if (on == T_FALSE)
            status = GR_DOMAIN;
        else if (on == T_UNKNOWN)
            status = GR_UNABLE;
    }

    if (status != GR_SUCCESS)
        res->is_infinity = T_UNKNOWN;

    return status;
}

int
gr_ec_aff_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return GR_DOMAIN;
    if (P->is_infinity == T_UNKNOWN)
        return GR_UNABLE;

    status |= gr_set(x, AX(P), R);
    status |= gr_set(y, AY(P), R);

    return status;
}

int
gr_ec_aff_point_lift_x(gr_ec_aff_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t, u, y;
    int status = GR_SUCCESS;

    GR_TMP_INIT3(t, u, y, R);

    if (gr_ec_ctx_model(ctx) == GR_EC_SHORT_WEIERSTRASS)
    {
        /* y^2 = x^3 + a4 x + a6 */
        status |= gr_sqr(t, x, R);
        status |= gr_mul(t, t, x, R);
        status |= gr_mul(u, GR_EC_A4(ctx), x, R);
        status |= gr_add(t, t, u, R);
        status |= gr_add(t, t, GR_EC_A6(ctx), R);

        if (status == GR_SUCCESS)
            status = gr_sqrt(y, t, R);
    }
    else
    {
        /* (2y + a1 x + a3)^2 = 4 x^3 + b2 x^2 + 2 b4 x + b6 */
        status |= gr_sqr(t, x, R);
        status |= gr_mul(u, t, x, R);
        status |= gr_mul_ui(u, u, 4, R);
        status |= gr_mul(t, t, GR_EC_B2(ctx), R);
        status |= gr_add(u, u, t, R);
        status |= gr_mul(t, GR_EC_B4(ctx), x, R);
        status |= gr_mul_two(t, t, R);
        status |= gr_add(u, u, t, R);
        status |= gr_add(u, u, GR_EC_B6(ctx), R);

        if (status == GR_SUCCESS)
            status = gr_sqrt(y, u, R);

        /* y = (w - a1 x - a3) / 2 */
        if (status == GR_SUCCESS)
        {
            status |= gr_mul(t, GR_EC_A1(ctx), x, R);
            status |= gr_sub(y, y, t, R);
            status |= gr_sub(y, y, GR_EC_A3(ctx), R);
            status |= gr_set_si(t, 2, R);
            if (status == GR_SUCCESS)
                status = gr_div(y, y, t, R);
        }
    }

    if (status == GR_SUCCESS)
        status = _gr_ec_aff_point_set_affine(res, x, y, ctx);
    else
        res->is_infinity = T_UNKNOWN;

    GR_TMP_CLEAR3(t, u, y, R);

    return status;
}

truth_t
gr_ec_aff_point_equal(const gr_ec_aff_point_t P, const gr_ec_aff_point_t Q,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);

    if (P->is_infinity == T_UNKNOWN || Q->is_infinity == T_UNKNOWN)
        return T_UNKNOWN;

    if (P->is_infinity == T_TRUE || Q->is_infinity == T_TRUE)
        return (P->is_infinity == Q->is_infinity) ? T_TRUE : T_FALSE;

    return truth_and(gr_equal(AX(P), AX(Q), R), gr_equal(AY(P), AY(Q), R));
}

truth_t
gr_ec_aff_point_is_on_curve(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t res;
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return T_TRUE;
    if (P->is_infinity == T_UNKNOWN)
        return T_UNKNOWN;

    GR_TMP_INIT_VEC(t, 3, R);

    /* T0 = y^2 + a1 x y + a3 y */
    status |= gr_sqr(T(0), AY(P), R);
    status |= gr_mul(T(1), AX(P), AY(P), R);
    status |= gr_addmul(T(0), T(1), GR_EC_A1(ctx), R);
    status |= gr_addmul(T(0), AY(P), GR_EC_A3(ctx), R);

    /* T0 -= x^3 + a2 x^2 + a4 x + a6 */
    status |= gr_sqr(T(1), AX(P), R);
    status |= gr_mul(T(2), T(1), AX(P), R);
    status |= gr_sub(T(0), T(0), T(2), R);
    status |= gr_submul(T(0), T(1), GR_EC_A2(ctx), R);
    status |= gr_submul(T(0), AX(P), GR_EC_A4(ctx), R);
    status |= gr_sub(T(0), T(0), GR_EC_A6(ctx), R);

    res = (status == GR_SUCCESS) ? gr_is_zero(T(0), R) : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 3, R);

    return res;
}

int
gr_ec_aff_point_write(gr_stream_t out, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return gr_stream_write(out, "O");
    if (P->is_infinity == T_UNKNOWN)
        return gr_stream_write(out, "(invalid point)");

    status |= gr_stream_write(out, "(");
    status |= gr_write(out, AX(P), R);
    status |= gr_stream_write(out, ", ");
    status |= gr_write(out, AY(P), R);
    status |= gr_stream_write(out, ")");

    return status;
}

int
gr_ec_aff_point_print(const gr_ec_aff_point_t P, gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_ec_aff_point_write(out, P, ctx);
}

int
gr_ec_aff_point_randtest(gr_ec_aff_point_t res, flint_rand_t state,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr x, y;
    int status = GR_UNABLE;
    slong iter;

    if (n_randint(state, 8) == 0)
        return gr_ec_aff_point_zero(res, ctx);

    GR_TMP_INIT2(x, y, R);

    for (iter = 0; iter < 20; iter++)
    {
        status = gr_randtest(x, state, R);

        if (status == GR_SUCCESS)
            status = gr_ec_aff_point_lift_x(res, x, ctx);

        if (status != GR_DOMAIN)
            break;
    }

    /* lifting divides by 2 in residue characteristic 2, so over a finite
       ring fall back to trying random pairs */
    if (status != GR_SUCCESS && gr_ctx_is_finite(R) == T_TRUE)
    {
        for (iter = 0; iter < 100; iter++)
        {
            if (gr_randtest(x, state, R) != GR_SUCCESS
                    || gr_randtest(y, state, R) != GR_SUCCESS)
                break;

            if (_gr_ec_aff_point_set_affine(res, x, y, ctx) == GR_SUCCESS
                    && gr_ec_aff_point_is_on_curve(res, ctx) == T_TRUE)
            {
                status = GR_SUCCESS;
                break;
            }
        }
    }

    if (status != GR_SUCCESS)
    {
        /* failing to find a point does not mean there is none */
        if (status == GR_DOMAIN)
            status = GR_UNABLE;

        res->is_infinity = T_UNKNOWN;
    }

    GR_TMP_CLEAR2(x, y, R);

    return status;
}
