/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_ec.h"

#define PX(P) GR_EC_POINT_X(P, ctx)
#define PY(P) GR_EC_POINT_Y(P, ctx)
#define PZ(P) GR_EC_POINT_Z(P, ctx)

#define AX(P) GR_EC_AFF_POINT_X(P, ctx)
#define AY(P) GR_EC_AFF_POINT_Y(P, ctx)

#define JX(P) GR_EC_JAC_POINT_X(P, ctx)
#define JY(P) GR_EC_JAC_POINT_Y(P, ctx)
#define JZ(P) GR_EC_JAC_POINT_Z(P, ctx)

int
gr_ec_point_set_aff_point(gr_ec_point_t res, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    if (P->is_infinity == T_TRUE)
        return gr_ec_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
        return GR_UNABLE;

    return _gr_ec_point_set_affine(res, AX(P), AY(P), ctx);
}

int
gr_ec_point_set_jac_point(gr_ec_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (P->is_infinity == T_TRUE)
        return gr_ec_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
        return GR_UNABLE;

    /* (X : Y : Z) with x = X/Z^2, y = Y/Z^3 is (X Z : Y : Z^3) */
    status |= gr_sqr(PZ(res), JZ(P), R);
    status |= gr_mul(PZ(res), PZ(res), JZ(P), R);
    status |= gr_mul(PX(res), JX(P), JZ(P), R);
    status |= gr_set(PY(res), JY(P), R);

    return status;
}

int
gr_ec_jac_point_set_point(gr_ec_jac_point_t res, const gr_ec_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    truth_t inf = gr_ec_point_is_inf(P, ctx);
    int status = GR_SUCCESS;

    if (inf == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);
    if (inf == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    /* affine (X/Z, Y/Z) is the Jacobian triple (X Z, Y Z^2, Z) */
    status |= gr_sqr(JY(res), PZ(P), R);
    status |= gr_mul(JY(res), JY(res), PY(P), R);
    status |= gr_mul(JX(res), PX(P), PZ(P), R);
    status |= gr_set(JZ(res), PZ(P), R);

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}

int
gr_ec_jac_point_set_aff_point(gr_ec_jac_point_t res, const gr_ec_aff_point_t P,
        gr_ec_ctx_t ctx)
{
    if (P->is_infinity == T_TRUE)
        return gr_ec_jac_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    return _gr_ec_jac_point_set_affine(res, AX(P), AY(P), ctx);
}

int
gr_ec_aff_point_set_point(gr_ec_aff_point_t res, const gr_ec_point_t P,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    truth_t inf = gr_ec_point_is_inf(P, ctx);
    gr_ptr t;
    int status = GR_SUCCESS;

    if (inf == T_TRUE)
        return gr_ec_aff_point_zero(res, ctx);
    if (inf == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    GR_TMP_INIT(t, R);

    status = gr_inv(t, PZ(P), R);

    if (status == GR_SUCCESS)
    {
        status |= gr_mul(AX(res), PX(P), t, R);
        status |= gr_mul(AY(res), PY(P), t, R);
    }

    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    GR_TMP_CLEAR(t, R);

    return status;
}

int
gr_ec_aff_point_set_jac_point(gr_ec_aff_point_t res, const gr_ec_jac_point_t P,
        gr_ec_ctx_t ctx)
{
    int status;

    if (P->is_infinity == T_TRUE)
        return gr_ec_aff_point_zero(res, ctx);
    if (P->is_infinity == T_UNKNOWN)
    {
        res->is_infinity = T_UNKNOWN;
        return GR_UNABLE;
    }

    status = gr_ec_jac_point_get_affine(AX(res), AY(res), P, ctx);
    res->is_infinity = (status == GR_SUCCESS) ? T_FALSE : T_UNKNOWN;

    return status;
}
