/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void gcd_check(
    fq_nmod_mpoly_t g,
    fq_nmod_mpoly_t a,
    fq_nmod_mpoly_t b,
    fq_nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fq_nmod_mpoly_t ca, cb, cg;

    fq_nmod_mpoly_init(ca, ctx);
    fq_nmod_mpoly_init(cb, ctx);
    fq_nmod_mpoly_init(cg, ctx);

    res = fq_nmod_mpoly_gcd_brown(g, a, b, ctx);
    fq_nmod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (fq_nmod_mpoly_is_zero(g, ctx))
    {
        if (!fq_nmod_mpoly_is_zero(a, ctx) || !fq_nmod_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fq_nmod_is_one(g->coeffs + 0, ctx->fqctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd is monic\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = 1;
    res = res && fq_nmod_mpoly_divides(ca, a, g, ctx);
    res = res && fq_nmod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = fq_nmod_mpoly_gcd_brown(cg, ca, cb, ctx);
    fq_nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (!fq_nmod_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

cleanup:

    fq_nmod_mpoly_clear(ca, ctx);
    fq_nmod_mpoly_clear(cb, ctx);
    fq_nmod_mpoly_clear(cg, ctx);
}

int
main(void)
{
    slong tmul = 5;
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....");
    fflush(stdout);

    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t g, a, b;
        const char * vars[] = {"x", "y", "z"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 3, ORD_LEX, 2, 2);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_set_str_pretty(a, "(x+y+z^2)*(x-y^9+z^3)", vars, ctx);
        fq_nmod_mpoly_set_str_pretty(b, "(x+y+z^9)*(x^9+y+z^2)", vars, ctx);

        gcd_check(g, a, b, ctx, 0, 0, "example");

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul*flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g;
        slong len, len1, len2;
        slong degbound;
        flint_bitcnt_t pbits;
        slong deg;

        pbits = 1 + n_randint(state, FLINT_BITS);
        pbits = 1 + n_randint(state, pbits);
        deg = 1 + n_randint(state, 4);
        fq_nmod_mpoly_ctx_init_rand(ctx, state, 4, pbits, deg);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 150);
        len2 = n_randint(state, 150);

        degbound = 1 + 50/ctx->minfo->nvars/ctx->minfo->nvars;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(g, state, len, degbound, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fq_nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fq_nmod_mpoly_mul(a, a, g, ctx);
            fq_nmod_mpoly_mul(b, b, g, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, ctx, i, j, "random dense");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
