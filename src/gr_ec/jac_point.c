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

#define JX(P) GR_EC_JAC_POINT_X(P, ctx)
#define JY(P) GR_EC_JAC_POINT_Y(P, ctx)
#define JZ(P) GR_EC_JAC_POINT_Z(P, ctx)

#define T(i) GR_ENTRY(t, i, sz)

void
gr_ec_jac_point_init(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);

    P->coords = flint_malloc(3 * R->sizeof_elem);
    _gr_vec_init(P->coords, 3, R);
    P->is_infinity = T_TRUE;
}

void
gr_ec_jac_point_clear(gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
{
    _gr_vec_clear(P->coords, 3, GR_EC_ELEM_CTX(ctx));
    flint_free(P->coords);
    P->coords = NULL;
}

int
gr_ec_jac_point_set(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (P->is_infinity != T_TRUE)
        status = _gr_vec_set(res->coords, P->coords, 3, GR_EC_ELEM_CTX(ctx));

    res->is_infinity = P->is_infinity;

    return (P->is_infinity == T_UNKNOWN) ? GR_UNABLE : status;
}

int
gr_ec_jac_point_zero(gr_ec_jac_point_t res, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    res->is_infinity = T_TRUE;
    return GR_SUCCESS;
}

int
_gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(JX(res), x, R);
    status |= gr_set(JY(res), y, R);
    status |= gr_one(JZ(res), R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

int
gr_ec_jac_point_set_affine(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y,
        gr_ec_ctx_t ctx)
{
    int status = _gr_ec_jac_point_set_affine(res, x, y, ctx);

    if (status == GR_SUCCESS)
    {
        truth_t on = gr_ec_jac_point_is_on_curve(res, ctx);

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
_gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y,
        gr_srcptr z, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(JX(res), x, R);
    status |= gr_set(JY(res), y, R);
    status |= gr_set(JZ(res), z, R);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

int
gr_ec_jac_point_set_jacobian(gr_ec_jac_point_t res, gr_srcptr x, gr_srcptr y,
        gr_srcptr z, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = _gr_ec_jac_point_set_jacobian(res, x, y, z, ctx);

    if (status == GR_SUCCESS && gr_is_zero(z, R) == T_TRUE)
        status = GR_DOMAIN;

    if (status == GR_SUCCESS)
    {
        truth_t on = gr_ec_jac_point_is_on_curve(res, ctx);

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
gr_ec_jac_point_get_affine(gr_ptr x, gr_ptr y, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr zi, zz;
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return GR_DOMAIN;
    if (P->is_infinity == T_UNKNOWN)
        return GR_UNABLE;

    GR_TMP_INIT2(zi, zz, R);

    status |= gr_inv(zi, JZ(P), R);

    if (status == GR_SUCCESS)
    {
        status |= gr_sqr(zz, zi, R);
        status |= gr_mul(x, JX(P), zz, R);
        status |= gr_mul(zz, zz, zi, R);
        status |= gr_mul(y, JY(P), zz, R);
    }

    GR_TMP_CLEAR2(zi, zz, R);

    return status;
}

int
gr_ec_jac_point_normalize(gr_ec_jac_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ptr x, y;
    int status;

    if (P->is_infinity == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
        return GR_UNABLE;

    GR_TMP_INIT2(x, y, R);

    status = gr_ec_jac_point_get_affine(x, y, P, ctx);

    if (status == GR_SUCCESS)
        status = _gr_ec_jac_point_set_affine(res, x, y, ctx);
    else
        res->is_infinity = T_UNKNOWN;

    GR_TMP_CLEAR2(x, y, R);

    return status;
}

truth_t
gr_ec_jac_point_equal(const gr_ec_jac_point_t P, const gr_ec_jac_point_t Q,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    slong sz = R->sizeof_elem;
    gr_ptr t;
    truth_t res;
    int status = GR_SUCCESS;

    if (P->is_infinity == T_UNKNOWN || Q->is_infinity == T_UNKNOWN)
        return T_UNKNOWN;

    if (P->is_infinity == T_TRUE || Q->is_infinity == T_TRUE)
        return (P->is_infinity == Q->is_infinity) ? T_TRUE : T_FALSE;

    GR_TMP_INIT_VEC(t, 4, R);

    /* X1 Z2^2 == X2 Z1^2 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_sqr(T(1), JZ(Q), R);
    status |= gr_mul(T(2), JX(P), T(1), R);
    status |= gr_mul(T(3), JX(Q), T(0), R);
    status |= gr_sub(T(2), T(2), T(3), R);
    res = gr_is_zero(T(2), R);

    if (res != T_FALSE)
    {
        /* Y1 Z2^3 == Y2 Z1^3 */
        status |= gr_mul(T(0), T(0), JZ(P), R);
        status |= gr_mul(T(1), T(1), JZ(Q), R);
        status |= gr_mul(T(2), JY(P), T(1), R);
        status |= gr_mul(T(3), JY(Q), T(0), R);
        status |= gr_sub(T(2), T(2), T(3), R);
        res = truth_and(res, gr_is_zero(T(2), R));
    }

    GR_TMP_CLEAR_VEC(t, 4, R);

    return (status == GR_SUCCESS) ? res : T_UNKNOWN;
}

truth_t
gr_ec_jac_point_is_on_curve(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
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

    GR_TMP_INIT_VEC(t, 6, R);

    /* T0 = Z^2, T1 = Z^3, T2 = Z^4, T3 = Z^6 */
    status |= gr_sqr(T(0), JZ(P), R);
    status |= gr_mul(T(1), T(0), JZ(P), R);
    status |= gr_sqr(T(2), T(0), R);
    status |= gr_mul(T(3), T(2), T(0), R);

    /* T4 = Y^2 + a1 X Y Z + a3 Y Z^3 */
    status |= gr_sqr(T(4), JY(P), R);
    status |= gr_mul(T(5), JX(P), JY(P), R);
    status |= gr_mul(T(5), T(5), JZ(P), R);
    status |= gr_addmul(T(4), T(5), GR_EC_A1(ctx), R);
    status |= gr_mul(T(5), JY(P), T(1), R);
    status |= gr_addmul(T(4), T(5), GR_EC_A3(ctx), R);

    /* T4 -= X^3 + a2 X^2 Z^2 + a4 X Z^4 + a6 Z^6 */
    status |= gr_sqr(T(5), JX(P), R);
    status |= gr_mul(T(1), T(5), JX(P), R);
    status |= gr_sub(T(4), T(4), T(1), R);
    status |= gr_mul(T(5), T(5), T(0), R);
    status |= gr_submul(T(4), T(5), GR_EC_A2(ctx), R);
    status |= gr_mul(T(5), JX(P), T(2), R);
    status |= gr_submul(T(4), T(5), GR_EC_A4(ctx), R);
    status |= gr_submul(T(4), T(3), GR_EC_A6(ctx), R);

    res = (status == GR_SUCCESS) ? gr_is_zero(T(4), R) : T_UNKNOWN;

    GR_TMP_CLEAR_VEC(t, 6, R);

    return res;
}

int
gr_ec_jac_point_write(gr_stream_t out, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return gr_stream_write(out, "O");
    if (P->is_infinity == T_UNKNOWN)
        return gr_stream_write(out, "(invalid point)");

    status |= gr_stream_write(out, "(");
    status |= gr_write(out, JX(P), R);
    status |= gr_stream_write(out, ", ");
    status |= gr_write(out, JY(P), R);
    status |= gr_stream_write(out, ", ");
    status |= gr_write(out, JZ(P), R);
    status |= gr_stream_write(out, ")");

    return status;
}

int
gr_ec_jac_point_get_str(char ** res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    int status;
    gr_stream_init_str(out);
    status = gr_ec_jac_point_write(out, P, ctx);
    *res = out->s;
    return status;
}

int
gr_ec_jac_point_print(const gr_ec_jac_point_t P, gr_ec_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_ec_jac_point_write(out, P, ctx);
}

int
gr_ec_jac_point_randtest(gr_ec_jac_point_t res, flint_rand_t state,
        gr_ec_ctx_t ctx)
{
    gr_ec_point_t P;
    int status;

    gr_ec_point_init(P, ctx);

    status = gr_ec_point_randtest(P, state, ctx);

    if (status == GR_SUCCESS)
        status = gr_ec_jac_point_set_point(res, P, ctx);

    if (status != GR_SUCCESS)
        res->is_infinity = T_UNKNOWN;

    gr_ec_point_clear(P, ctx);

    return status;
}
