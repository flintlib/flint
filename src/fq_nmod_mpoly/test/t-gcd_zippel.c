/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

/* Defined in t-gcd_brown.c, t-gcd.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_zippel2.c, t-gcd_zippel.c */
#define gcd_check gcd_check_gcd_zippel
void gcd_check(
    fq_nmod_mpoly_t g,
    fq_nmod_mpoly_t a,
    fq_nmod_mpoly_t b,
    fq_nmod_mpoly_t t,
    fq_nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    fq_nmod_mpoly_t ca, cb, cg;

    fq_nmod_mpoly_init(ca, ctx);
    fq_nmod_mpoly_init(cb, ctx);
    fq_nmod_mpoly_init(cg, ctx);

    if (!fq_nmod_mpoly_gcd_zippel(g, a, b, ctx))
    {
        flint_printf("FAIL: check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fq_nmod_mpoly_assert_canonical(g, ctx);

    if (fq_nmod_mpoly_is_zero(g, ctx))
    {
        if (!fq_nmod_mpoly_is_zero(a, ctx) || !fq_nmod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fq_nmod_mpoly_is_monic(g, ctx))
    {
        flint_printf("FAIL: check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fq_nmod_mpoly_is_zero(t, ctx) && !fq_nmod_mpoly_divides(cg, g, t, ctx))
    {
        flint_printf("FAIL: check gcd divisor\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fq_nmod_mpoly_divides(ca, a, g, ctx) ||
        !fq_nmod_mpoly_divides(cb, b, g, ctx))
    {
        flint_printf("FAIL: check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fq_nmod_mpoly_gcd_zippel(cg, ca, cb, ctx))
    {
        flint_printf("FAIL: check cofactor gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fq_nmod_mpoly_assert_canonical(cg, ctx);

    if (!fq_nmod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fq_nmod_mpoly_clear(ca, ctx);
    fq_nmod_mpoly_clear(cb, ctx);
    fq_nmod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fq_nmod_mpoly_gcd_zippel, state)
{
    slong i, j;

    for (i = 0; i < 20*flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds, * degbounds1, * degbounds2;
        flint_bitcnt_t pbits;
        slong deg;

        pbits = 1 + n_randint(state, FLINT_BITS);
        pbits = 1 + n_randint(state, pbits);
        deg = 2 + n_randint(state, 4);
        fq_nmod_mpoly_ctx_init_rand(ctx, state, 5, pbits, deg);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(ca, ctx);
        fq_nmod_mpoly_init(cb, ctx);
        fq_nmod_mpoly_init(cg, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len = n_randint(state, 30) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 50/(2*ctx->minfo->nvars - 1);
        degbounds = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds1 = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds2 = (ulong * ) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds[j] = n_randint(state, degbound + 1) + 1;
            degbounds1[j] = n_randint(state, degbound + 1) + 1;
            degbounds2[j] = n_randint(state, degbound + 1) + 1;
        }

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bounds(t, state, len, degbounds, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bounds(a, state, len1, degbounds1, ctx);
            fq_nmod_mpoly_randtest_bounds(b, state, len2, degbounds2, ctx);

            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "random");
        }

        flint_free(degbounds);
        flint_free(degbounds1);
        flint_free(degbounds2);

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(ca, ctx);
        fq_nmod_mpoly_clear(cb, ctx);
        fq_nmod_mpoly_clear(cg, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
