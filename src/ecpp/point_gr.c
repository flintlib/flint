/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
#include "gr_ec.h"
#include "mpn_mod.h"
#include "ecpp.h"

/*
    Point arithmetic on y^2 = x^3 + a x + b over Z/n, in Jacobian
    coordinates with Z = 0 for the point at infinity.

    The arithmetic itself is gr_ec's. What this module needs beyond the
    ordinary group law is the guarantee that the law was entitled to the
    branches it took: the formulas are polynomial and always compute
    something, but the choice between adding, doubling and the point at
    infinity is made by testing whether a quantity vanishes, and modulo a
    composite n a quantity can vanish for one prime factor and not
    another. The branch is then right for one factor and wrong for the
    other, and the result is quietly meaningless -- which for a primality
    proof is the one thing that must not happen.

    gr_ec_jac_point_mul_fmpz_witness collects those quantities into acc,
    so that gcd(acc, n) != 1 says the run is not to be trusted; see the
    notes there. That guarantee is the whole of what this file adds.

    Only a4 is handed over because the Jacobian formulas never read a6.
*/

#define PX(P) (P)
#define PY(P) GR_ENTRY(P, 1, sz)
#define PZ(P) GR_ENTRY(P, 2, sz)

void
ecpp_gr_ctx_init(gr_ctx_t gctx, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);

    if (fmpz_size(n) >= MPN_MOD_MIN_LIMBS && fmpz_size(n) <= MPN_MOD_MAX_LIMBS
            && gr_ctx_init_mpn_mod(gctx, n) == GR_SUCCESS)
        return;

    gr_ctx_init_fmpz_mod(gctx, n);
}

int
ecpp_point_mul_gr(gr_ptr R, gr_srcptr P, const fmpz_t k, gr_srcptr a,
        gr_ptr acc, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    gr_ec_ctx_t E;
    gr_ec_jac_point_t gP, gR;
    gr_ptr zero;
    int ok;

    zero = gr_heap_init(ctx);

    /* a6 plays no part in the group law, so any value will do */
    ok = (gr_zero(zero, ctx) == GR_SUCCESS)
        && (gr_ec_ctx_init_short_weierstrass(E, ctx, a, zero) == GR_SUCCESS);

    gr_heap_clear(zero, ctx);

    if (!ok)
        return 0;

    gr_ec_jac_point_init(gP, E);
    gr_ec_jac_point_init(gR, E);

    ok = (gr_set(GR_EC_JAC_POINT_X(gP, E), PX(P), ctx) == GR_SUCCESS)
        && (gr_set(GR_EC_JAC_POINT_Y(gP, E), PY(P), ctx) == GR_SUCCESS)
        && (gr_set(GR_EC_JAC_POINT_Z(gP, E), PZ(P), ctx) == GR_SUCCESS);

    if (ok)
    {
        gP->is_infinity = gr_is_zero(PZ(P), ctx);
        ok = (gP->is_infinity != T_UNKNOWN);
    }

    if (ok)
        ok = (gr_ec_jac_point_mul_fmpz_witness(gR, acc, gP, k, E) == GR_SUCCESS);

    if (ok)
    {
        if (gr_ec_jac_point_is_inf(gR, E) == T_TRUE)
            ok = (gr_zero(PZ(R), ctx) == GR_SUCCESS);
        else
            ok = (gr_set(PX(R), GR_EC_JAC_POINT_X(gR, E), ctx) == GR_SUCCESS)
                && (gr_set(PY(R), GR_EC_JAC_POINT_Y(gR, E), ctx) == GR_SUCCESS)
                && (gr_set(PZ(R), GR_EC_JAC_POINT_Z(gR, E), ctx) == GR_SUCCESS);
    }

    gr_ec_jac_point_clear(gR, E);
    gr_ec_jac_point_clear(gP, E);
    gr_ec_ctx_clear(E);

    return ok && gr_is_zero(acc, ctx) != T_TRUE;
}

/* the fmpz_mod interface of point.c, through gr */
int
ecpp_point_mul(ecpp_point_t R, const ecpp_point_t P, const fmpz_t k,
                        const fmpz_t a, fmpz_t acc, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gctx;
    gr_ptr gP, gR, ga, gacc;
    slong sz;
    int ok;

    ecpp_gr_ctx_init(gctx, ctx);
    sz = gctx->sizeof_elem;
    gP = gr_heap_init_vec(3, gctx);
    gR = gr_heap_init_vec(3, gctx);
    ga = gr_heap_init(gctx);
    gacc = gr_heap_init(gctx);

    ok = (gr_set_fmpz(PX(gP), P->X, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(PY(gP), P->Y, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(PZ(gP), P->Z, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(ga, a, gctx) == GR_SUCCESS)
      && (gr_set_fmpz(gacc, acc, gctx) == GR_SUCCESS);

    if (ok)
        ok = ecpp_point_mul_gr(gR, gP, k, ga, gacc, gctx);

    /* acc and R are read back even on failure (acc may be zero) */
    GR_MUST_SUCCEED(gr_get_fmpz(acc, gacc, gctx));

    if (ok)
    {
        GR_MUST_SUCCEED(gr_get_fmpz(R->X, PX(gR), gctx));
        GR_MUST_SUCCEED(gr_get_fmpz(R->Y, PY(gR), gctx));
        GR_MUST_SUCCEED(gr_get_fmpz(R->Z, PZ(gR), gctx));
    }
    else
        fmpz_zero(R->Z);

    gr_heap_clear_vec(gP, 3, gctx);
    gr_heap_clear_vec(gR, 3, gctx);
    gr_heap_clear(ga, gctx);
    gr_heap_clear(gacc, gctx);
    gr_ctx_clear(gctx);

    return ok;
}
