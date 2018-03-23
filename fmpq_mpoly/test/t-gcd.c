/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_t lc;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound;
        slong coeff_bits;
        int res;

        fmpq_init(lc);

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(ca, ctx);
        fmpq_mpoly_init(cb, ctx);
        fmpq_mpoly_init(cg, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 20/(1 + ctx->zctx->minfo->nvars);

        coeff_bits = n_randint(state, 20);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (fmpq_mpoly_is_zero(t, ctx));
            fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);

            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bound(g, state, len, coeff_bits, degbound, ctx);
            res = fmpq_mpoly_gcd(g, a, b, ctx);
            if (!res) {
                continue;
            }

            fmpq_mpoly_assert_canonical(g, ctx);

            if (fmpq_mpoly_is_zero(g, ctx))
            {
                if (!fmpq_mpoly_is_zero(a, ctx) || !fmpq_mpoly_is_zero(b, ctx)) {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            fmpq_mpoly_get_coeff_fmpq(lc, g, 0, ctx);
            if (!fmpq_is_one(lc))
            {
                printf("FAIL\n");
                flint_printf("Check gcd is monic\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && fmpq_mpoly_divides(ca, a, g, ctx);
            res = res && fmpq_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
            fmpq_mpoly_assert_canonical(ca, ctx);
            fmpq_mpoly_assert_canonical(cb, ctx);

            fmpq_mpoly_gcd(cg, ca, cb, ctx);
            fmpq_mpoly_assert_canonical(cg, ctx);

            if (!fmpq_mpoly_equal_fmpq(cg, lc, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(ca, ctx);
        fmpq_mpoly_clear(cb, ctx);
        fmpq_mpoly_clear(cg, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);

        fmpq_clear(lc);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

