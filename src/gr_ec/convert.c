/*
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Opus 5 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_ec.h"

#define AX(P) GR_EC_AFF_POINT_X(P, ctx)
#define AY(P) GR_EC_AFF_POINT_Y(P, ctx)

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
