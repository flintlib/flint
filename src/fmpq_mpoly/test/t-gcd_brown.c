/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

/* Defined in t-gcd.c, t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_subresultant.c, t-gcd_zippel2.c */
#define gcd_check gcd_check_gcd_brown
void gcd_check(
    fmpq_mpoly_t g,
    fmpq_mpoly_t a,
    fmpq_mpoly_t b,
    fmpq_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fmpq_mpoly_t ca, cb, cg;

    fmpq_mpoly_init(ca, ctx);
    fmpq_mpoly_init(cb, ctx);
    fmpq_mpoly_init(cg, ctx);

    res = fmpq_mpoly_gcd_brown(g, a, b, ctx);
    fmpq_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (fmpq_mpoly_is_zero(g, ctx))
    {
        if (!fmpq_mpoly_is_zero(a, ctx) || !fmpq_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpq_mpoly_is_monic(g, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd has positive lc\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && fmpq_mpoly_divides(ca, a, g, ctx);
    res = res && fmpq_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpq_mpoly_gcd_brown(cg, ca, cb, ctx);
    fmpq_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpq_mpoly_clear(ca, ctx);
    fmpq_mpoly_clear(cb, ctx);
    fmpq_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fmpq_mpoly_gcd_brown, state)
{
    slong tmul = 10;
    slong i, j;

    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t g, a, b;
        const char * vars[] = {"x", "y", "z"};

        fmpq_mpoly_ctx_init(ctx, 3, ORD_LEX);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_set_str_pretty(a, "(x+y+z^2)*(x-y^9+z^3)", vars, ctx);
        fmpq_mpoly_set_str_pretty(b, "(x+y+z^2)*(x^9+y+z^2)", vars, ctx);

        gcd_check(g, a, b, ctx, 0, 0, "example");

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);

        len = n_randint(state, 40) + 1;
        len1 = n_randint(state, 60);
        len2 = n_randint(state, 60);

        n = FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx));
        degbound = 1 + 50/n/n;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(g, state, len, coeff_bits + 1, degbound, ctx);
            if (fmpq_mpoly_is_zero(g, ctx))
                fmpq_mpoly_one(g, ctx);
            fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpq_mpoly_mul(a, a, g, ctx);
            fmpq_mpoly_mul(b, b, g, ctx);
            fmpq_mpoly_scalar_mul_ui(a, a, n_randint(state, 10) + 1, ctx);
            fmpq_mpoly_scalar_mul_ui(b, b, n_randint(state, 10) + 1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            gcd_check(g, a, b, ctx, i, j, "random dense");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
