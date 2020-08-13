/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

void gcd_check(
    fmpz_mpoly_t g,
    fmpz_mpoly_t a,
    fmpz_mpoly_t b,
    fmpz_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fmpz_mpoly_t ca, cb, cg;

    fmpz_mpoly_init(ca, ctx);
    fmpz_mpoly_init(cb, ctx);
    fmpz_mpoly_init(cg, ctx);

    res = fmpz_mpoly_gcd_brown_threaded(g, a, b, ctx);
    fmpz_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (fmpz_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            flint_abort();
        }
        goto cleanup;
    }

    if (fmpz_sgn(g->coeffs + 0) <= 0)
    {
        printf("FAIL\n");
        flint_printf("Check gcd has positive lc\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mpoly_divides(ca, a, g, ctx);
    res = res && fmpz_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = fmpz_mpoly_gcd_brown_threaded(cg, ca, cb, ctx);
    fmpz_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (!fmpz_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

cleanup:

    fmpz_mpoly_clear(ca, ctx);
    fmpz_mpoly_clear(cb, ctx);
    fmpz_mpoly_clear(cg, ctx);
}


int
main(void)
{
    slong i, j;
    slong tmul = 10;
    slong max_threads = 5;
    FLINT_TEST_INIT(state);
#ifdef _WIN32
    tmul = 1;
#endif

    flint_printf("gcd_brown_threaded....");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b;
        const char * vars[] = {"x", "y", "z"};

        fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);

        fmpz_mpoly_set_str_pretty(a, "x+123412341*y+z^2", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "123497345963*x^3+y+z", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "x+y^4+3459687345*z", vars, ctx);
        fmpz_mpoly_mul(a, a, g, ctx);
        fmpz_mpoly_mul(b, b, g, ctx);
        fmpz_mpoly_zero(g, ctx);

        gcd_check(g, a, b, ctx, -1, -1, "example");

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);

        len = n_randint(state, 40) + 1;
        len1 = n_randint(state, 80);
        len2 = n_randint(state, 80);

        degbound = 1 + 50/ctx->minfo->nvars/ctx->minfo->nvars;

        coeff_bits = n_randint(state, 500);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(g, state, len, coeff_bits + 1, degbound, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, g, ctx);
            fmpz_mpoly_mul(b, b, g, ctx);
            fmpz_mpoly_scalar_mul_ui(a, a, n_randint(state, 10) + 1, ctx);
            fmpz_mpoly_scalar_mul_ui(b, b, n_randint(state, 10) + 1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            gcd_check(g, a, b, ctx, i, j, "random dense");
        }

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

