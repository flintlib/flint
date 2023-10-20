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

TEST_FUNCTION_START(fq_nmod_mpoly_derivative, state)
{
    int i, j, result;
    slong tmul = 5;

    /* Check d(f*g) = df*g + f*dg */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, fp, gp, hp, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        slong idx;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        if (fq_nmod_mpoly_ctx_nvars(ctx) < 1)
        {
            fq_nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(fp, ctx);
        fq_nmod_mpoly_init(gp, ctx);
        fq_nmod_mpoly_init(hp, ctx);
        fq_nmod_mpoly_init(t1, ctx);
        fq_nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        fq_nmod_mpoly_randtest_bits(hp, state, len, exp_bits, ctx);
        fq_nmod_mpoly_randtest_bits(fp, state, len, exp_bits1, ctx);
        fq_nmod_mpoly_randtest_bits(gp, state, len, exp_bits2, ctx);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);

            idx = n_randint(state, fq_nmod_mpoly_ctx_nvars(ctx));

            fq_nmod_mpoly_mul(h, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            fq_nmod_mpoly_derivative(hp, h, idx, ctx);
            fq_nmod_mpoly_assert_canonical(hp, ctx);

            fq_nmod_mpoly_derivative(fp, f, idx, ctx);
            fq_nmod_mpoly_assert_canonical(fp, ctx);
            fq_nmod_mpoly_derivative(gp, g, idx, ctx);
            fq_nmod_mpoly_assert_canonical(gp, ctx);

            fq_nmod_mpoly_mul(t1, f, gp, ctx);
            fq_nmod_mpoly_assert_canonical(t1, ctx);
            fq_nmod_mpoly_mul(t2, g, fp, ctx);
            fq_nmod_mpoly_assert_canonical(t2, ctx);
            fq_nmod_mpoly_add(t1, t1, t2, ctx);
            fq_nmod_mpoly_assert_canonical(t1, ctx);

            result = fq_nmod_mpoly_equal(hp, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg\n"
                                                  "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(fp, ctx);
        fq_nmod_mpoly_clear(gp, ctx);
        fq_nmod_mpoly_clear(hp, ctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);

        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check d(f*g) = df*g + f*dg with aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, fp, gp, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        slong idx;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 2);
        if (fq_nmod_mpoly_ctx_nvars(ctx) < 1)
        {
            fq_nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(fp, ctx);
        fq_nmod_mpoly_init(gp, ctx);
        fq_nmod_mpoly_init(t1, ctx);
        fq_nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
        fq_nmod_mpoly_randtest_bits(fp, state, len, exp_bits, ctx);
        fq_nmod_mpoly_randtest_bits(gp, state, len, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            idx = n_randint(state, fq_nmod_mpoly_ctx_nvars(ctx));

            fq_nmod_mpoly_mul(h, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            fq_nmod_mpoly_derivative(h, h, idx, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_set(fp, f, ctx);
            fq_nmod_mpoly_derivative(fp, fp, idx, ctx);
            fq_nmod_mpoly_assert_canonical(fp, ctx);
            fq_nmod_mpoly_set(gp, g, ctx);
            fq_nmod_mpoly_derivative(gp, gp, idx, ctx);
            fq_nmod_mpoly_assert_canonical(gp, ctx);

            fq_nmod_mpoly_mul(t1, f, gp, ctx);
            fq_nmod_mpoly_assert_canonical(t1, ctx);
            fq_nmod_mpoly_mul(t2, g, fp, ctx);
            fq_nmod_mpoly_assert_canonical(t2, ctx);
            fq_nmod_mpoly_add(t1, t1, t2, ctx);
            fq_nmod_mpoly_assert_canonical(t1, ctx);

            result = fq_nmod_mpoly_equal(h, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg with aliasing\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(fp, ctx);
        fq_nmod_mpoly_clear(gp, ctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);

        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
