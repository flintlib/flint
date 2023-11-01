/*
    Copyright (C) 2019 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly_factor.h"

/* Defined in t-gcd_brown.c, t-gcd_brown_threaded.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#define compute_gcd compute_gcd_zippel2
int compute_gcd(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, max_deg;
    flint_bitcnt_t wbits;
    int success = 0;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mpoly_t Ac, Bc, Gc, Gamma, lcAl, lcBl;
    slong * Adegs, * Bdegs, * perm;
    ulong * shift, * stride;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
            fmpz_mpoly_zero(G, ctx);
        else if (fmpz_sgn(B->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, B, ctx);
        else
            fmpz_mpoly_set(G, B, ctx);
        return 1;
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, A, ctx);
        else
            fmpz_mpoly_set(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mpoly_gcd_zippel(G, A, B, ctx);
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(ctx->minfo->nvars >= 3);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    Adegs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    Bdegs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    perm = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    shift = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);
    stride = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);

    mpoly_degrees_si(Adegs, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(Bdegs, B->exps, B->length, B->bits, ctx->minfo);

    max_deg = 0;
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
        max_deg = FLINT_MAX(max_deg, Adegs[i]);
        max_deg = FLINT_MAX(max_deg, Bdegs[i]);
    }

    fmpz_mpoly_ctx_init(lctx, ctx->minfo->nvars, ORD_LEX);

    wbits = 1 + FLINT_BIT_COUNT(max_deg);
    wbits = FLINT_MAX(MPOLY_MIN_BITS, wbits);
    wbits = mpoly_fix_bits(wbits, lctx->minfo);
    FLINT_ASSERT(wbits <= FLINT_BITS);

    fmpz_mpoly_init3(Al, A->length, wbits, lctx);
    fmpz_mpoly_init3(Bl, B->length, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gamma, 0, wbits, lctx);
    fmpz_mpoly_init3(lcAl, 0, wbits, lctx);
    fmpz_mpoly_init3(lcBl, 0, wbits, lctx);

    if (FLINT_BIT_COUNT(FLINT_MAX(Adegs[0], Adegs[1])) >= FLINT_BITS/2 ||
        FLINT_BIT_COUNT(FLINT_MAX(Bdegs[0], Bdegs[1])) >= FLINT_BITS/2)
    {
        success = 0;
        goto cleanup;
    }

    fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);
    fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, perm, shift, stride);

    success = fmpz_mpolyl_content(Ac, Al, 2, lctx) &&
              fmpz_mpolyl_content(Bc, Bl, 2, lctx);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_gcd(Gc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_divides(Al, Al, Ac, lctx);
    FLINT_ASSERT(success);

    success = fmpz_mpoly_divides(Bl, Bl, Bc, lctx);
    FLINT_ASSERT(success);

    fmpz_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    fmpz_mpolyl_lead_coeff(lcAl, Al, 2, lctx);
    fmpz_mpolyl_lead_coeff(lcBl, Bl, 2, lctx);
    success = fmpz_mpoly_gcd(Gamma, lcAl, lcBl, lctx);
    if (!success)
        goto cleanup;

    success = fmpz_mpolyl_gcd_zippel2(Gl, Abarl, Bbarl, Al, Bl, Gamma, lctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mpoly_from_mpolyl_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                                Gl, lctx, perm, shift, stride);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    success = 1;

cleanup:

    flint_free(Adegs);
    flint_free(Bdegs);
    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_clear(Ac, lctx);
    fmpz_mpoly_clear(Bc, lctx);
    fmpz_mpoly_clear(Gc, lctx);
    fmpz_mpoly_clear(Gamma, lctx);
    fmpz_mpoly_clear(lcAl, lctx);
    fmpz_mpoly_clear(lcBl, lctx);

    fmpz_mpoly_ctx_clear(lctx);

    return success;
}

/* Defined in t-gcd_brown.c, t-gcd_brown_threaded.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#ifndef gcd_check
#define gcd_check gcd_check
void gcd_check(
    fmpz_mpoly_t g,
    fmpz_mpoly_t a,
    fmpz_mpoly_t b,
    const fmpz_mpoly_t gdiv,
    fmpz_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name,
    int compute_gcd_fun(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t))
{
    int res;
    fmpz_mpoly_t ca, cb, cg;

    fmpz_mpoly_init(ca, ctx);
    fmpz_mpoly_init(cb, ctx);
    fmpz_mpoly_init(cg, ctx);

    res = compute_gcd_fun(g, a, b, ctx);

    fmpz_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpz_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (fmpz_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (fmpz_sgn(g->coeffs + 0) <= 0)
    {
        flint_printf("FAIL: Check gcd has positive lc\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mpoly_divides(ca, a, g, ctx);
    res = res && fmpz_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = compute_gcd_fun(cg, ca, cb, ctx);
    fmpz_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpz_mpoly_clear(ca, ctx);
    fmpz_mpoly_clear(cb, ctx);
    fmpz_mpoly_clear(cg, ctx);
}
#endif

TEST_FUNCTION_START(fmpz_mpoly_factor_gcd_zippel2, state)
{
    slong i, j, tmul = 20;

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b, t;
        const char* vars[] = {"y", "t", "x", "z"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        fmpz_mpoly_set_str_pretty(t, "x+y+z+t", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "x^2+y^2+z^2+t^2", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "x^3+y^3+z^3+t^3", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 0, "example", compute_gcd);

        fmpz_mpoly_set_str_pretty(t, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 0, "example", compute_gcd);

        fmpz_mpoly_set_str_pretty(t, "y + t^2 + x^3 + z^4", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "y*t + 1 + (x - z^5)*(y + t)", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "y*t + 1 + (x - z^5)*(y - t + x)", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 1, "trigger unlucky ksub", compute_gcd);

        fmpz_mpoly_set_str_pretty(t, "y + 33857*t^2 + 35153*x^3 + 40433*z^4", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "y^4 + t^3 + x^2 + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "y^3 + t^4 + x^2 + z", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 2, "trigger zipple no match", compute_gcd);

        fmpz_mpoly_set_str_pretty(t, "y + t + x^3 + z^3", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "(x - z^4)*y + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(x + z)*y + t", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 3, "trigger ksub lc kill", compute_gcd);

        fmpz_mpoly_set_str_pretty(t, "y + t + x^3 + z^3", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "(x - z^4 + 33857*(x*z))*y + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(x + z)*y + t", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 4, "trigger ksub lc kill mod p", compute_gcd);

        fmpz_mpoly_set_str_pretty(t, "(1 + x^10)*t*y + t + x + z", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "(1 + x + t)*(x - 33857*x^2 + 35153*z^4)*y*t + z + x*y", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(2*x - t^2)*(x - 33857*x^2 + 35153*z^4)*y^2*t + t*z + y", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 5, "trigger gcd lc terms vanish mod p", compute_gcd);

        if (FLINT_BITS == 64)
        {
            fmpz_mpoly_set_str_pretty(t, "t*y + t + x^9999999999 + z^9999999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(a, "t + y + x^9999999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "t^2 + t*z^9999999999 + y + 1", vars, ctx);
        }
        else
        {
            fmpz_mpoly_set_str_pretty(t, "t*y + t + x^999999 + z^9999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(a, "t + y + x^999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "t^2 + t*z^999999 + y + 1", vars, ctx);
        }
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 6, "trigger big p", compute_gcd);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        if (ctx->minfo->nvars < 3)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 35) + 1;
        len1 = n_randint(state, 35) + 1;
        len2 = n_randint(state, 35) + 1;

        degbound = 2 + 100/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            if (fmpz_mpoly_is_zero(t, ctx))
                fmpz_mpoly_one(t, ctx);

            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "sparse", compute_gcd);
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef compute_gcd
