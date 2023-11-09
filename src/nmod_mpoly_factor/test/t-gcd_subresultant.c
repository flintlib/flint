/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly_factor.h"

/* Defined in t-gcd_subresultant.c and t-gcd_zippel.c */
#define compute_gcd compute_gcd_subresultant
int compute_gcd(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong var = 0;
    nmod_mpoly_t Ac, Bc, Gc, s, t;
    nmod_mpoly_univar_t Ax, Bx, Gx;

    if (nmod_mpoly_is_zero(A, ctx) || nmod_mpoly_is_zero(B, ctx))
        return nmod_mpoly_gcd(G, A, B, ctx);

    nmod_mpoly_init(Ac, ctx);
    nmod_mpoly_init(Bc, ctx);
    nmod_mpoly_init(Gc, ctx);
    nmod_mpoly_init(s, ctx);
    nmod_mpoly_init(t, ctx);
    nmod_mpoly_univar_init(Ax, ctx);
    nmod_mpoly_univar_init(Bx, ctx);
    nmod_mpoly_univar_init(Gx, ctx);

    nmod_mpoly_to_univar(Ax, A, var, ctx);
    nmod_mpoly_to_univar(Bx, B, var, ctx);

    success = _nmod_mpoly_vec_content_mpoly(Ac, Ax->coeffs, Ax->length, ctx) &&
              _nmod_mpoly_vec_content_mpoly(Bc, Bx->coeffs, Bx->length, ctx) &&
              nmod_mpoly_gcd(Gc, Ac, Bc, ctx);
    if (!success)
        goto cleanup;

    _nmod_mpoly_vec_divexact_mpoly(Ax->coeffs, Ax->length, Ac, ctx);
    _nmod_mpoly_vec_divexact_mpoly(Bx->coeffs, Bx->length, Bc, ctx);

    success = fmpz_cmp(Ax->exps + 0, Bx->exps + 0) > 0 ?
                            nmod_mpoly_univar_pseudo_gcd(Gx, Ax, Bx, ctx) :
                            nmod_mpoly_univar_pseudo_gcd(Gx, Bx, Ax, ctx);
    if (!success)
        goto cleanup;

    if (nmod_mpoly_gcd(t, Ax->coeffs + 0, Bx->coeffs + 0, ctx) &&
                                                                t->length == 1)
    {
        nmod_mpoly_term_content(s, Gx->coeffs + 0, ctx);
        nmod_mpoly_divexact(t, Gx->coeffs + 0, s, ctx);
        _nmod_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    }
    else if (nmod_mpoly_gcd(t, Ax->coeffs + Ax->length - 1,
                           Bx->coeffs + Bx->length - 1, ctx) && t->length == 1)
    {
        nmod_mpoly_term_content(s, Gx->coeffs + Gx->length - 1, ctx);
        nmod_mpoly_divexact(t, Gx->coeffs + Gx->length - 1, s, ctx);
        _nmod_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    }

    success = _nmod_mpoly_vec_content_mpoly(t, Gx->coeffs, Gx->length, ctx);
    if (!success)
        goto cleanup;

    _nmod_mpoly_vec_divexact_mpoly(Gx->coeffs, Gx->length, t, ctx);
    _nmod_mpoly_vec_mul_mpoly(Gx->coeffs, Gx->length, Gc, ctx);
    _nmod_mpoly_from_univar(G, FLINT_MIN(A->bits, B->bits), Gx, var, ctx);
    nmod_mpoly_make_monic(G, G, ctx);

    success = 1;

cleanup:

    nmod_mpoly_clear(Ac, ctx);
    nmod_mpoly_clear(Bc, ctx);
    nmod_mpoly_clear(Gc, ctx);
    nmod_mpoly_clear(s, ctx);
    nmod_mpoly_clear(t, ctx);
    nmod_mpoly_univar_clear(Ax, ctx);
    nmod_mpoly_univar_clear(Bx, ctx);
    nmod_mpoly_univar_clear(Gx, ctx);

    return success;
}

/* Defined in t-gcd_subresultant.c and t-gcd_zippel.c */
#define gcd_check gcd_check_gcd_subresultant
void gcd_check(
    nmod_mpoly_t g,
    nmod_mpoly_t a,
    nmod_mpoly_t b,
    const nmod_mpoly_t gdiv,
    nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    nmod_mpoly_t ca, cb, cg;

    nmod_mpoly_init(ca, ctx);
    nmod_mpoly_init(cb, ctx);
    nmod_mpoly_init(cg, ctx);

    res = compute_gcd(g, a, b, ctx);

    nmod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_is_zero(gdiv, ctx))
    {
        if (!nmod_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (nmod_mpoly_is_zero(g, ctx))
    {
        if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (g->coeffs[0] != 1)
    {
        flint_printf("FAIL: Check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && nmod_mpoly_divides(ca, a, g, ctx);
    res = res && nmod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = compute_gcd(cg, ca, cb, ctx);
    nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    nmod_mpoly_clear(ca, ctx);
    nmod_mpoly_clear(cb, ctx);
    nmod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(nmod_mpoly_factor_gcd_subresultant, state)
{
    slong i, j, tmul = 15;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);
        if (ctx->minfo->nvars < 1)
        {
            nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 12) + 1;
        len1 = n_randint(state, 12);
        len2 = n_randint(state, 12);

        degbound = 1 + 10/ctx->minfo->nvars;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);
            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "random dense");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef compute_gcd
#undef gcd_check
