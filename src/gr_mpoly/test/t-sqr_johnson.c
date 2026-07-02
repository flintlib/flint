/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2026 sqr_johnson author

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_sqr_johnson, state)
{
    slong i, j;

    /* Check sqr_johnson(f) == mul_johnson(f, f) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, k1, k2;
        slong len1;
        flint_bitcnt_t exp_bits1;
        int status;

        gr_ctx_init_random(cctx, state);

        /* gr_mpoly_sqr_johnson assumes a (approximately) commutative
           coefficient ring a priori; skip rings that are not known to be so */
        if (gr_ctx_is_commutative_ring(cctx) != T_TRUE &&
            gr_ctx_is_approx_commutative_ring(cctx) != T_TRUE)
        {
            gr_ctx_clear(cctx);
            continue;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(k1, ctx);
        gr_mpoly_init(k2, ctx);

        if (gr_ctx_is_finite(cctx) == T_TRUE)
            len1 = n_randint(state, 100);
        else
            len1 = n_randint(state, 5);

        exp_bits1 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(k1, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(k2, state, len1, exp_bits1, ctx);

            /* k1 = f^2 via dedicated squaring */
            status |= gr_mpoly_sqr_johnson(k1, f, ctx);

            if (status == GR_SUCCESS)
                gr_mpoly_assert_canonical(k1, ctx);

            /* k2 = f*f via generic multiplication (reference) */
            status |= gr_mpoly_mul_johnson(k2, f, f, ctx);

            if (status == GR_SUCCESS)
                gr_mpoly_assert_canonical(k2, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(k1, k2, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check sqr(f) == f * f\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("sqr(f)  = "); gr_mpoly_print_pretty(k1, ctx); flint_printf("\n");
                flint_printf("f * f   = "); gr_mpoly_print_pretty(k2, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* Also exercise aliasing: f = f^2 in place */
        {
            status = GR_SUCCESS;
            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);

            status |= gr_mpoly_set(k2, f, ctx);
            status |= gr_mpoly_mul_johnson(k2, k2, k2, ctx);

            status |= gr_mpoly_sqr_johnson(f, f, ctx);

            if (status == GR_SUCCESS)
                gr_mpoly_assert_canonical(f, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(f, k2, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check aliased sqr(f) == f * f\n");
                flint_printf("i = %wd\n", i);
                gr_ctx_println(ctx);
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(k1, ctx);
        gr_mpoly_clear(k2, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
