/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Multiplying a point by an element of another ring.

    The group of points is a Z-module and nothing more, so an element of an
    arbitrary ring has no business acting on it. There is one exception
    worth making: if m kills every point then k P depends only on k modulo
    m, so Z/m acts on the curve, and Z/m is exactly the ring a caller
    reaches for when the scalars are known modulo the order of the group.

    The condition is that m annihilate the group, which the context is
    asked about through gr_ec_ctx_annihilates. Notably this does not count
    the points: if the order is already known the answer is a divisibility
    test, and if it is not, the modulus is tested against random points.
    Handing over a ring of scalars should not cost a point count.

    Reading the scalar out needs gr_get_fmpz in the scalar ring, which
    nmod, fmpz_mod and mpn_mod all provide.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/* the modulus of a Z/n, or GR_DOMAIN if y_ctx is not one */
static int
_scalar_ring_modulus(fmpz_t n, gr_ctx_t y_ctx)
{
    switch (y_ctx->which_ring)
    {
        case GR_CTX_NMOD:
        case GR_CTX_NMOD8:
        case GR_CTX_NMOD32:
        case GR_CTX_FMPZ_MOD:
        case GR_CTX_MPN_MOD:
            return gr_ctx_cardinality_fmpz(n, y_ctx);

        default:
            return GR_DOMAIN;
    }
}

/*
    The integer a foreign scalar stands for, once its ring has been checked
    against the curve. GR_DOMAIN if the ring cannot act on this curve,
    GR_UNABLE if that could not be decided.
*/
int
_gr_ec_scalar_of_other(fmpz_t k, gr_srcptr y, gr_ctx_t y_ctx, gr_ec_ctx_t ctx)
{
    fmpz_t n;
    truth_t acts;
    int status;

    fmpz_init(n);

    status = _scalar_ring_modulus(n, y_ctx);

    if (status == GR_SUCCESS)
    {
        acts = gr_ec_ctx_annihilates(ctx, n);

        if (acts == T_FALSE)
            status = GR_DOMAIN;
        else if (acts == T_UNKNOWN)
            status = GR_UNABLE;
    }

    if (status == GR_SUCCESS)
        status = gr_get_fmpz(k, y, y_ctx);

    fmpz_clear(n);

    return status;
}

#define GR_EC_MUL_OTHER_IMPL(kind) \
int \
_ ## kind ## _mul_other(kind ## _t res, const kind ## _t P, gr_srcptr y, \
        gr_ctx_t y_ctx, gr_ec_ctx_t ctx) \
{ \
    fmpz_t k; \
    int status; \
 \
    fmpz_init(k); \
    status = _gr_ec_scalar_of_other(k, y, y_ctx, ctx); \
 \
    if (status == GR_SUCCESS) \
        status = kind ## _mul_fmpz(res, P, k, ctx); \
 \
    fmpz_clear(k); \
 \
    return status; \
} \
 \
int \
_ ## kind ## _other_mul(kind ## _t res, gr_srcptr y, gr_ctx_t y_ctx, \
        const kind ## _t P, gr_ec_ctx_t ctx) \
{ \
    return _ ## kind ## _mul_other(res, P, y, y_ctx, ctx); \
}

GR_EC_MUL_OTHER_IMPL(gr_ec_point)
GR_EC_MUL_OTHER_IMPL(gr_ec_aff_point)
GR_EC_MUL_OTHER_IMPL(gr_ec_jac_point)
