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

TEST_FUNCTION_START(fmpq_mpoly_gcd_zippel, state)
{
    slong i, j;
    slong tmul = 10;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, ca, cb, cg, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds;
        int res;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(ca, ctx);
        fmpq_mpoly_init(cb, ctx);
        fmpq_mpoly_init(cg, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 15) + 1;
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 15);

        degbound = 100/(2*fmpq_mpoly_ctx_nvars(ctx) - 1);
        degbounds = FLINT_ARRAY_ALLOC(fmpq_mpoly_ctx_nvars(ctx), ulong);
        for (j = 0; j < fmpq_mpoly_ctx_nvars(ctx); j++)
            degbounds[j] = n_randint(state, degbound + UWORD(1)) + UWORD(1);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bounds(t, state, len, coeff_bits + 1, degbounds, ctx);
            if (fmpq_mpoly_is_zero(g, ctx))
                fmpq_mpoly_one(g, ctx);
            fmpq_mpoly_randtest_bounds(a, state, len1, coeff_bits, degbounds, ctx);
            fmpq_mpoly_randtest_bounds(b, state, len2, coeff_bits, degbounds, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            res = fmpq_mpoly_gcd_zippel(g, a, b, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check that gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            if (fmpq_mpoly_is_zero(g, ctx))
            {
                if (!fmpq_mpoly_is_zero(a, ctx) || !fmpq_mpoly_is_zero(b, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }
                continue;
            }

            if (!fmpq_mpoly_is_monic(g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check gcd has positive lc\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            res = 1;
            res = res && fmpq_mpoly_divides(ca, a, g, ctx);
            res = res && fmpq_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            res = fmpq_mpoly_gcd_zippel(cg, ca, cb, ctx);

            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check that cofactor gcd could be computed\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            if (!fmpq_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(degbounds);

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(ca, ctx);
        fmpq_mpoly_clear(cb, ctx);
        fmpq_mpoly_clear(cg, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
