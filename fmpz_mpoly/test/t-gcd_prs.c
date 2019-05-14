/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"

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

    res = fmpz_mpoly_gcd_prs(g, a, b, ctx);
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

    res = fmpz_mpoly_gcd_prs(cg, ca, cb, ctx);
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
    int i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_prs....");
    fflush(stdout);

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g;
        ordering_t ord;
        slong nvars, len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 3) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);

        for (j = 0; j < 4; j++)
        {
            len1 = n_randint(state, 18);
            len2 = n_randint(state, 18);
            len3 = n_randint(state, 14);
            exp_bound1 = n_randint(state, 4) + 2;
            exp_bound2 = n_randint(state, 4) + 2;
            exp_bound3 = n_randint(state, 3) + 2;
            coeff_bits = n_randint(state, 3) + 1;
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(g, state, len3, coeff_bits, exp_bound3, ctx);

            fmpz_mpoly_mul_johnson(a, a, g, ctx);
            fmpz_mpoly_mul_johnson(b, b, g, ctx);

            fmpz_mpoly_randtest_bound(g, state, len3, coeff_bits, exp_bound3, ctx);

            gcd_check(g, a, b, ctx, i, j, "random small");
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

