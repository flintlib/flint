/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly_factor.h"

/* Defined in t-gcd_brown.c, t-gcd_hensel.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#define compute_gcd compute_gcd_brown
int compute_gcd(
    fmpz_mod_mpoly_t G,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong nvars = ctx->minfo->nvars;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    flint_bitcnt_t wbits;
    fmpz_mod_mpoly_ctx_t nctx;
    fmpz_mod_mpolyn_t An, Bn, Gn, Abarn, Bbarn;
    fmpz_mod_poly_polyun_mpolyn_stack_t St;

    if (fmpz_mod_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mod_mpoly_is_zero(B, ctx))
            fmpz_mod_mpoly_zero(G, ctx);
        else
            fmpz_mod_mpoly_make_monic(G, B, ctx);

        return 1;
    }

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (nvars < 2)
    {
        return fmpz_mod_mpoly_gcd(G, A, B, ctx);
    }

    perm = (slong *) flint_malloc(nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(nvars*sizeof(ulong));
    for (i = 0; i < nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mod_mpoly_ctx_init(nctx, nvars, ORD_LEX, fmpz_mod_mpoly_ctx_modulus(ctx));
    fmpz_mod_mpolyn_init(An, wbits, nctx);
    fmpz_mod_mpolyn_init(Bn, wbits, nctx);
    fmpz_mod_mpolyn_init(Gn, wbits, nctx);
    fmpz_mod_mpolyn_init(Abarn, wbits, nctx);
    fmpz_mod_mpolyn_init(Bbarn, wbits, nctx);
    fmpz_mod_poly_stack_init(St->poly_stack);
    fmpz_mod_polyun_stack_init(St->polyun_stack);
    fmpz_mod_mpolyn_stack_init(St->mpolyn_stack, wbits, nctx);

    fmpz_mod_mpoly_to_mpolyn_perm_deflate(An, nctx, A, ctx, perm, shift, stride);
    fmpz_mod_mpoly_to_mpolyn_perm_deflate(Bn, nctx, B, ctx, perm, shift, stride);

    success = fmpz_mod_mpolyn_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                    nvars - 1, nctx, NULL, St);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_from_mpolyn_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                                Gn, nctx, perm, shift, stride);
    fmpz_mod_mpoly_make_monic(G, G, ctx);

    success = 1;

cleanup:

    fmpz_mod_poly_stack_clear(St->poly_stack);
    fmpz_mod_polyun_stack_clear(St->polyun_stack);
    fmpz_mod_mpolyn_stack_clear(St->mpolyn_stack, nctx);
    fmpz_mod_mpolyn_clear(An, nctx);
    fmpz_mod_mpolyn_clear(Bn, nctx);
    fmpz_mod_mpolyn_clear(Gn, nctx);
    fmpz_mod_mpolyn_clear(Abarn, nctx);
    fmpz_mod_mpolyn_clear(Bbarn, nctx);
    fmpz_mod_mpoly_ctx_clear(nctx);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}

/* Defined in t-gcd_brown.c, t-gcd_hensel.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#ifndef gcd_check
#define gcd_check gcd_check
void gcd_check(
    fmpz_mod_mpoly_t g,
    fmpz_mod_mpoly_t a,
    fmpz_mod_mpoly_t b,
    const fmpz_mod_mpoly_t gdiv,
    fmpz_mod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name,
    int (* compute_gcd_fun)(fmpz_mod_mpoly_t, const fmpz_mod_mpoly_t, const fmpz_mod_mpoly_t, const fmpz_mod_mpoly_ctx_t))
{
    int res;
    fmpz_mod_mpoly_t ca, cb, cg;

    fmpz_mod_mpoly_init(ca, ctx);
    fmpz_mod_mpoly_init(cb, ctx);
    fmpz_mod_mpoly_init(cg, ctx);

    res = compute_gcd_fun(g, a, b, ctx);

    fmpz_mod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        if (fmpz_bits(fmpz_mod_mpoly_ctx_modulus(ctx)) < 10)
            goto cleanup;

        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpz_mod_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (fmpz_mod_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(a, ctx) || !fmpz_mod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpz_is_one(g->coeffs + 0))
    {
        flint_printf("FAIL: Check gcd has positive lc\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mod_mpoly_divides(ca, a, g, ctx);
    res = res && fmpz_mod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = compute_gcd_fun(cg, ca, cb, ctx);
    fmpz_mod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        if (fmpz_bits(fmpz_mod_mpoly_ctx_modulus(ctx)) < 10)
            goto cleanup;

        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpz_mod_mpoly_clear(ca, ctx);
    fmpz_mod_mpoly_clear(cb, ctx);
    fmpz_mod_mpoly_clear(cg, ctx);
}
#endif

TEST_FUNCTION_START(fmpz_mod_mpoly_factor_gcd_brown, state)
{
    slong i, j, tmul = 5;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong n, degbound;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 5, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len = n_randint(state, 40) + 1;
        len1 = n_randint(state, 60);
        len2 = n_randint(state, 60);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        degbound = 1 + 50/n/n;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "random dense", compute_gcd);
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef compute_gcd
