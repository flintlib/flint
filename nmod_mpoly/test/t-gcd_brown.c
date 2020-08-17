/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void gcd_check(
    nmod_mpoly_t g,
    nmod_mpoly_t a,
    nmod_mpoly_t b,
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

    nmod_mpoly_assert_canonical(a, ctx);
    nmod_mpoly_assert_canonical(b, ctx);
    res = nmod_mpoly_gcd_brown(g, a, b, ctx);
    nmod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (nmod_mpoly_is_zero(g, ctx))
    {
        if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            flint_abort();
        }
        goto cleanup;
    }

    if (g->coeffs[0] != UWORD(1))
    {
        printf("FAIL\n");
        flint_printf("Check gcd is monic\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = 1;
    res = res && nmod_mpoly_divides(ca, a, g, ctx);
    res = res && nmod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    nmod_mpoly_assert_canonical(ca, ctx);
    nmod_mpoly_assert_canonical(cb, ctx);
    res = nmod_mpoly_gcd_brown(cg, ca, cb, ctx);
    nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (!nmod_mpoly_equal_ui(cg, 1, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

cleanup:

    nmod_mpoly_clear(ca, ctx);
    nmod_mpoly_clear(cb, ctx);
    nmod_mpoly_clear(cg, ctx);
}


int
main(void)
{
    slong tmul = 5;
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_brown....");
    fflush(stdout);

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t p;

        p = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        p = n_randbits(state, p);
        p = n_nextprime(p, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, p < 3000 ? 4 : 5, p);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        degbound = 1 + 80/ctx->minfo->nvars/ctx->minfo->nvars;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(g, state, len, degbound, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul(a, a, g, ctx);
            nmod_mpoly_mul(b, b, g, ctx);
            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, ctx, i, j, "random dense");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

