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
#include "gr_vec.h"
#include "gr_ec.h"

#define PX(P) GR_EC_POINT_X(P, ctx)
#define PY(P) GR_EC_POINT_Y(P, ctx)
#define PZ(P) GR_EC_POINT_Z(P, ctx)

void
gr_ec_point_init(gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);

    P->coords = flint_malloc(3 * R->sizeof_elem);
    _gr_vec_init(P->coords, 3, R);
    GR_MUST_SUCCEED(gr_one(PY(P), R));
}

void
gr_ec_point_clear(gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    _gr_vec_clear(P->coords, 3, GR_EC_ELEM_CTX(ctx));
    flint_free(P->coords);
    P->coords = NULL;
}

int
gr_ec_point_set(gr_ec_point_t res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    return _gr_vec_set(res->coords, P->coords, 3, GR_EC_ELEM_CTX(ctx));
}

int
gr_ec_point_zero(gr_ec_point_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_zero(PX(res), R);
    status |= gr_one(PY(res), R);
    status |= gr_zero(PZ(res), R);

    return status;
}

int
_gr_ec_point_set_affine(gr_ec_point_t res, gr_srcptr x, gr_srcptr y,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(PX(res), x, R);
    status |= gr_set(PY(res), y, R);
    status |= gr_one(PZ(res), R);

    return status;
}

int
gr_ec_point_set_affine(gr_ec_point_t res, gr_srcptr x, gr_srcptr y,
        gr_ec_ctx_t ctx)
{
    int status = _gr_ec_point_set_affine(res, x, y, ctx);

    if (status == GR_SUCCESS)
    {
        truth_t on = gr_ec_point_is_on_curve(res, ctx);

        if (on == T_FALSE)
            status = GR_DOMAIN;
        else if (on == T_UNKNOWN)
            status = GR_UNABLE;
    }

    return status;
}

int
_gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y,
        gr_srcptr z, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(PX(res), x, R);
    status |= gr_set(PY(res), y, R);
    status |= gr_set(PZ(res), z, R);

    return status;
}

int
gr_ec_point_set_projective(gr_ec_point_t res, gr_srcptr x, gr_srcptr y,
        gr_srcptr z, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = _gr_ec_point_set_projective(res, x, y, z, ctx);

    if (status == GR_SUCCESS)
    {
        truth_t zero = truth_and(truth_and(gr_is_zero(x, R), gr_is_zero(y, R)),
                                    gr_is_zero(z, R));

        if (zero == T_TRUE)
            status = GR_DOMAIN;
    }

    if (status == GR_SUCCESS)
    {
        truth_t on = gr_ec_point_is_on_curve(res, ctx);

        if (on == T_FALSE)
            status = GR_DOMAIN;
        else if (on == T_UNKNOWN)
            status = GR_UNABLE;
    }

    return status;
}

int
gr_ec_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    truth_t inf = gr_ec_point_is_inf(P, ctx);
    gr_ptr t;
    int status = GR_SUCCESS;

    if (inf == T_TRUE)
        return GR_DOMAIN;
    if (inf == T_UNKNOWN)
        return GR_UNABLE;

    GR_TMP_INIT(t, R);

    status |= gr_inv(t, PZ(P), R);

    if (status == GR_SUCCESS)
    {
        status |= gr_mul(x, PX(P), t, R);
        status |= gr_mul(y, PY(P), t, R);
    }

    GR_TMP_CLEAR(t, R);

    return status;
}

int
gr_ec_point_normalize(gr_ec_point_t res, const gr_ec_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    truth_t inf = gr_ec_point_is_inf(P, ctx);
    gr_ptr t;
    int status = GR_SUCCESS;

    if (inf == T_TRUE)
        return gr_ec_point_zero(res, ctx);
    if (inf == T_UNKNOWN)
        return GR_UNABLE;

    GR_TMP_INIT(t, R);

    status |= gr_inv(t, PZ(P), R);

    if (status == GR_SUCCESS)
    {
        status |= gr_mul(PX(res), PX(P), t, R);
        status |= gr_mul(PY(res), PY(P), t, R);
        status |= gr_one(PZ(res), R);
    }

    GR_TMP_CLEAR(t, R);

    return status;
}

int
gr_ec_point_lift_x(gr_ec_point_t res, gr_srcptr x, gr_ec_ctx_t ctx)
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
        status = _gr_ec_point_set_affine(res, x, y, ctx);

    GR_TMP_CLEAR3(t, u, y, R);

    return status;
}

truth_t
gr_ec_point_is_inf(const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    return gr_is_zero(PZ(P), GR_EC_ELEM_CTX(ctx));
}

truth_t
gr_ec_point_equal(const gr_ec_point_t P, const gr_ec_point_t Q,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr t, u;
    truth_t res;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(t, u, R);

    status |= gr_mul(t, PX(P), PY(Q), R);
    status |= gr_mul(u, PX(Q), PY(P), R);
    status |= gr_sub(t, t, u, R);
    res = gr_is_zero(t, R);

    if (res != T_FALSE)
    {
        status |= gr_mul(t, PY(P), PZ(Q), R);
        status |= gr_mul(u, PY(Q), PZ(P), R);
        status |= gr_sub(t, t, u, R);
        res = truth_and(res, gr_is_zero(t, R));
    }

    if (res != T_FALSE)
    {
        status |= gr_mul(t, PX(P), PZ(Q), R);
        status |= gr_mul(u, PX(Q), PZ(P), R);
        status |= gr_sub(t, t, u, R);
        res = truth_and(res, gr_is_zero(t, R));
    }

    GR_TMP_CLEAR2(t, u, R);

    return (status == GR_SUCCESS) ? res : T_UNKNOWN;
}

truth_t
gr_ec_point_is_on_curve(const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr acc, t, u, zz;
    truth_t res;
    int status = GR_SUCCESS;

    GR_TMP_INIT4(acc, t, u, zz, R);

    /* acc = Y^2 Z + a1 X Y Z + a3 Y Z^2 */
    status |= gr_sqr(acc, PY(P), R);
    status |= gr_mul(acc, acc, PZ(P), R);

    status |= gr_mul(t, PX(P), PY(P), R);
    status |= gr_mul(t, t, PZ(P), R);
    status |= gr_mul(t, t, GR_EC_A1(ctx), R);
    status |= gr_add(acc, acc, t, R);

    status |= gr_sqr(zz, PZ(P), R);
    status |= gr_mul(t, PY(P), zz, R);
    status |= gr_mul(t, t, GR_EC_A3(ctx), R);
    status |= gr_add(acc, acc, t, R);

    /* acc -= X^3 + a2 X^2 Z + a4 X Z^2 + a6 Z^3 */
    status |= gr_sqr(u, PX(P), R);
    status |= gr_mul(t, u, PX(P), R);
    status |= gr_sub(acc, acc, t, R);

    status |= gr_mul(t, u, PZ(P), R);
    status |= gr_mul(t, t, GR_EC_A2(ctx), R);
    status |= gr_sub(acc, acc, t, R);

    status |= gr_mul(t, PX(P), zz, R);
    status |= gr_mul(t, t, GR_EC_A4(ctx), R);
    status |= gr_sub(acc, acc, t, R);

    status |= gr_mul(t, zz, PZ(P), R);
    status |= gr_mul(t, t, GR_EC_A6(ctx), R);
    status |= gr_sub(acc, acc, t, R);

    res = (status == GR_SUCCESS) ? gr_is_zero(acc, R) : T_UNKNOWN;

    GR_TMP_CLEAR4(acc, t, u, zz, R);

    return res;
}

int
gr_ec_point_write(gr_stream_t out, const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_stream_write(out, "(");
    status |= gr_write(out, PX(P), R);
    status |= gr_stream_write(out, " : ");
    status |= gr_write(out, PY(P), R);
    status |= gr_stream_write(out, " : ");
    status |= gr_write(out, PZ(P), R);
    status |= gr_stream_write(out, ")");

    return status;
}

int
gr_ec_point_get_str(char ** res, const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    int status;
    gr_stream_init_str(out);
    status = gr_ec_point_write(out, P, ctx);
    *res = out->s;
    return status;
}

int
gr_ec_point_print(const gr_ec_point_t P, gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_ec_point_write(out, P, ctx);
}

int
gr_ec_point_randtest(gr_ec_point_t res, flint_rand_t state, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr x, y;
    int status = GR_UNABLE;
    slong iter;

    if (n_randint(state, 8) == 0)
        return gr_ec_point_zero(res, ctx);

    GR_TMP_INIT2(x, y, R);

    for (iter = 0; iter < 20; iter++)
    {
        status = gr_randtest(x, state, R);

        if (status == GR_SUCCESS)
            status = gr_ec_point_lift_x(res, x, ctx);

        if (status != GR_DOMAIN)
            break;
    }

    /*
        Over a small finite ring, fall back to trying random pairs. This is
        the only route in residue characteristic 2, where lifting an
        x-coordinate would have to divide by 2.
    */
    if (status != GR_SUCCESS && gr_ctx_is_finite(R) == T_TRUE)
    {
        for (iter = 0; iter < 100; iter++)
        {
            if (gr_randtest(x, state, R) != GR_SUCCESS
                    || gr_randtest(y, state, R) != GR_SUCCESS)
                break;

            if (_gr_ec_point_set_affine(res, x, y, ctx) == GR_SUCCESS
                    && gr_ec_point_is_on_curve(res, ctx) == T_TRUE)
            {
                status = GR_SUCCESS;
                break;
            }
        }
    }

    /* no point was found; that is a failure to produce one, not a
       statement that none exists */
    if (status == GR_DOMAIN)
        status = GR_UNABLE;

    /* randomize the projective representative */
    if (status == GR_SUCCESS)
    {
        if (gr_randtest(x, state, R) == GR_SUCCESS
                && gr_is_invertible(x, R) == T_TRUE)
        {
            status |= gr_mul(PX(res), PX(res), x, R);
            status |= gr_mul(PY(res), PY(res), x, R);
            status |= gr_mul(PZ(res), PZ(res), x, R);
        }
    }

    GR_TMP_CLEAR2(x, y, R);

    return status;
}
