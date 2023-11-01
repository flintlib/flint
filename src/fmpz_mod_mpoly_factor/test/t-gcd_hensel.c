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
#define compute_gcd compute_gcd_hensel
int compute_gcd(
    fmpz_mod_mpoly_t G,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    flint_rand_t state;
    int success = 0;
    fmpz_mod_mpoly_ctx_t lctx;
    fmpz_mod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mod_mpoly_t Ac, Bc, Gc;
    ulong * shift, * stride;
    slong * perm;
    slong Gdeg_bound;

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

    if (ctx->minfo->nvars < 3)
    {
        return fmpz_mod_mpoly_gcd(G, A, B, ctx);
    }

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(!fmpz_mod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mod_mpoly_is_zero(B, ctx));
    FLINT_ASSERT(ctx->minfo->nvars >= 3);

    flint_randinit(state);

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mod_mpoly_ctx_init(lctx, ctx->minfo->nvars, ORD_LEX,
                                              fmpz_mod_mpoly_ctx_modulus(ctx));

    fmpz_mod_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mod_mpoly_init3(Gc, 0, wbits, lctx);

    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);
    fmpz_mod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, perm, shift, stride);

    success = fmpz_mod_mpolyl_content(Ac, Al, 1, lctx) &&
              fmpz_mod_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    if (!fmpz_mod_mpoly_divides(Al, Al, Ac, lctx) ||
        !fmpz_mod_mpoly_divides(Bl, Bl, Bc, lctx))
    {
        flint_printf("FAIL: Check content divides\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mod_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    Gdeg_bound = FLINT_MIN(fmpz_mod_mpoly_degree_si(Al, 0, lctx),
                           fmpz_mod_mpoly_degree_si(Bl, 0, lctx));
    success = fmpz_mod_mpolyl_gcd_hensel_smprime(Gl, Gdeg_bound, Abarl, Bbarl,
                                                                 Al, Bl, lctx);
    if (!success)
        goto cleanup;

    success = fmpz_mod_mpoly_gcd(Gc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    fmpz_mod_mpoly_mul(Gl, Gl, Gc, lctx);

    fmpz_mod_mpoly_from_mpolyl_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                                Gl, lctx, perm, shift, stride);
    fmpz_mod_mpoly_make_monic(G, G, ctx);

    success = 1;

cleanup:

    fmpz_mod_mpoly_clear(Al, lctx);
    fmpz_mod_mpoly_clear(Bl, lctx);
    fmpz_mod_mpoly_clear(Gl, lctx);
    fmpz_mod_mpoly_clear(Abarl, lctx);
    fmpz_mod_mpoly_clear(Bbarl, lctx);
    fmpz_mod_mpoly_clear(Ac, lctx);
    fmpz_mod_mpoly_clear(Bc, lctx);
    fmpz_mod_mpoly_clear(Gc, lctx);

    fmpz_mod_mpoly_ctx_clear(lctx);

    flint_randclear(state);

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

TEST_FUNCTION_START(fmpz_mod_mpoly_factor_gcd_hensel, state)
{
    slong i, j, tmul = 10;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong degbound;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 5, (i & 1) ? 150 : 15);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len = n_randint(state, 30) + 1;
        len1 = n_randint(state, 30) + 1;
        len2 = n_randint(state, 30) + 1;

        degbound = 2 + 30/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);

            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "sparse", compute_gcd);
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
