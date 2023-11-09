/*
    Copyright (C) 2019 Daniel Schultz

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
#define gcd_check gcd_check_gcd_brown
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

    if (!fq_nmod_mpoly_gcd_brown(g, a, b, ctx))
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

    if (!fq_nmod_mpoly_gcd_brown(cg, ca, cb, ctx))
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

TEST_FUNCTION_START(fq_nmod_mpoly_gcd_brown, state)
{
    slong i, j;

    for (i = 0; i < 5*flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, t, g;
        slong len, len1, len2;
        slong n, degbound;
        flint_bitcnt_t pbits;
        slong deg;

        pbits = 1 + n_randint(state, FLINT_BITS);
        pbits = 1 + n_randint(state, pbits);
        deg = 1 + n_randint(state, 4);
        fq_nmod_mpoly_ctx_init_rand(ctx, state, 4, pbits, deg);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(t, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 150);
        len2 = n_randint(state, 150);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        degbound = 1 + 50/n/n;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fq_nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "random dense");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
