/*
    Copyright (C) 2026 Fredrik Johansson
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Fable 5.1 and Opus 5.5

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

void
ecpp_gr_ctx_init(gr_ctx_t gctx, const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);

    if (fmpz_size(n) >= MPN_MOD_MIN_LIMBS && fmpz_size(n) <= MPN_MOD_MAX_LIMBS
            && gr_ctx_init_mpn_mod(gctx, n) == GR_SUCCESS)
        return;

    gr_ctx_init_fmpz_mod(gctx, n);
}

/* 1 if gcd(x, n) = 1 */
static int
_is_unit(gr_srcptr x, const fmpz_t n, gr_ctx_t R)
{
    fmpz_t g;
    int res;

    fmpz_init(g);
    GR_MUST_SUCCEED(gr_get_fmpz(g, x, R));
    fmpz_gcd(g, g, n);
    res = fmpz_is_one(g);
    fmpz_clear(g);

    return res;
}

/*
    Modulo a composite n, the group law may take a wrong branch (a value
    can vanish modulo one prime factor of n but not another). The witness
    w collects every value treated as nonzero; the result is only trusted
    if w is a unit.
*/
int
ecpp_point_check(const fmpz_t x, const fmpz_t y, const fmpz_t a,
        const fmpz_t b, const fmpz_t k, const fmpz_t q,
        const fmpz_mod_ctx_t ctx)
{
    const fmpz * n = fmpz_mod_ctx_modulus(ctx);
    gr_ctx_t R;
    gr_ec_ctx_t E;
    gr_ec_jac_point_t P, Q;
    gr_ptr u, v, w;
    int result;

    ecpp_gr_ctx_init(R, ctx);
    u = gr_heap_init(R);
    v = gr_heap_init(R);
    w = gr_heap_init(R);

    GR_MUST_SUCCEED(gr_set_fmpz(u, a, R));
    GR_MUST_SUCCEED(gr_set_fmpz(v, b, R));

    /* a singular curve proves nothing */
    if (gr_ec_ctx_init_short_weierstrass(E, R, u, v) != GR_SUCCESS)
    {
        result = ECPP_POINT_WRONG_ORDER;
        goto cleanup_ring;
    }

    gr_ec_jac_point_init(P, E);
    gr_ec_jac_point_init(Q, E);

    GR_MUST_SUCCEED(gr_set_fmpz(u, x, R));
    GR_MUST_SUCCEED(gr_set_fmpz(v, y, R));
    GR_MUST_SUCCEED(_gr_ec_jac_point_set_affine(P, u, v, E));
    GR_MUST_SUCCEED(gr_one(w, R));

    /* Q = k P */
    if (gr_ec_jac_point_mul_fmpz_witness(Q, P, k, w, E) != GR_SUCCESS
            || gr_ec_jac_point_is_inf(Q, E) != T_FALSE)
    {
        result = ECPP_POINT_ZERO;
        goto cleanup;
    }

    if (!_is_unit(GR_EC_JAC_POINT_Z(Q, E), n, R))
    {
        result = ECPP_POINT_COMPOSITE;
        goto cleanup;
    }

    /* P = q Q */
    if (gr_ec_jac_point_mul_fmpz_witness(P, Q, q, w, E) != GR_SUCCESS)
        result = ECPP_POINT_WRONG_ORDER;
    else if (!_is_unit(w, n, R))
        result = ECPP_POINT_COMPOSITE;
    else if (gr_ec_jac_point_is_inf(P, E) != T_TRUE)
        result = ECPP_POINT_WRONG_ORDER;
    else
        result = ECPP_POINT_OK;

cleanup:
    gr_ec_jac_point_clear(P, E);
    gr_ec_jac_point_clear(Q, E);
    gr_ec_ctx_clear(E);

cleanup_ring:
    gr_heap_clear(u, R);
    gr_heap_clear(v, R);
    gr_heap_clear(w, R);
    gr_ctx_clear(R);

    return result;
}
