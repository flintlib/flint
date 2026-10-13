/*
    Copyright (C) 2026 Fredrik Johansson
    Copyright (C) 2026 Maël Hostettler
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
    Points in Jacobian coordinates (X : Y : Z), Z = 0 for the point at
    infinity.
*/

void
ecpp_point_init(ecpp_point_t P)
{
    fmpz_init(P->X);
    fmpz_init(P->Y);
    fmpz_init(P->Z);
}

void
ecpp_point_clear(ecpp_point_t P)
{
    fmpz_clear(P->X);
    fmpz_clear(P->Y);
    fmpz_clear(P->Z);
}

void
ecpp_point_set_affine(ecpp_point_t P, const fmpz_t x, const fmpz_t y)
{
    fmpz_set(P->X, x);
    fmpz_set(P->Y, y);
    fmpz_one(P->Z);
}

int
ecpp_point_is_zero(const ecpp_point_t P)
{
    return fmpz_is_zero(P->Z);
}

/*
    Point arithmetic on y^2 = x^3 + a x + b over Z/n, using gr_ec.

    Modulo a composite n, the group law may take a wrong branch (a value
    can be zero modulo one prime factor of n but not another), so the
    result is only trusted if every value treated as nonzero is a unit;
    gr_ec_jac_point_mul_fmpz_witness collects these values in acc.
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
    gr_ptr a6;
    int status;
    int ok;

    a6 = gr_heap_init(ctx);

    /*
        a6 plays no part in the Jacobian group law, but gr_ec rejects a
        singular curve: a6 = 0 is singular when a = 0, and then a6 = 1
        gives discriminant -432.
    */
    status = gr_zero(a6, ctx);
    if (status == GR_SUCCESS)
        status = gr_ec_ctx_init_short_weierstrass(E, ctx, a, a6);
    if (status == GR_DOMAIN && gr_one(a6, ctx) == GR_SUCCESS)
        status = gr_ec_ctx_init_short_weierstrass(E, ctx, a, a6);
    ok = (status == GR_SUCCESS);

    gr_heap_clear(a6, ctx);

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
        ok = (gr_ec_jac_point_mul_fmpz_witness(gR, gP, k, acc, E) == GR_SUCCESS);

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
