/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_mpoly_mul_monomial, state)
{
    slong i, j;

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h, k1, k2, t;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 4);

        gr_mpoly_init(f, mctx, cctx);
        gr_mpoly_init(g, mctx, cctx);
        gr_mpoly_init(h, mctx, cctx);
        gr_mpoly_init(k1, mctx, cctx);
        gr_mpoly_init(k2, mctx, cctx);
        gr_mpoly_init(t, mctx, cctx);

        if (cctx->methods == _ca_methods)
        {
            len = n_randint(state, 10);
            len1 = n_randint(state, 10);
        }
        else
        {
            len = n_randint(state, 100);
            len1 = n_randint(state, 100);
        }

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        /* todo: aliasing tests */

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, mctx, cctx);
            status |= gr_mpoly_randtest_bits(g, state, 1, exp_bits2, mctx, cctx);
            status |= gr_mpoly_randtest_bits(h, state, 1, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k1, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k2, state, len, exp_bits, mctx, cctx);

            if (g->length != 1 || h->length != 1)
                continue;

            status |= gr_mpoly_add(k1, g, h, mctx, cctx);

            if (n_randint(state, 2))
                status |= gr_mpoly_mul_johnson(k1, f, k1, mctx, cctx);
            else
                status |= gr_mpoly_mul_johnson(k1, k1, f, mctx, cctx);

            if (status == GR_SUCCESS)
                gr_mpoly_assert_canonical(k1, mctx, cctx);

            status |= gr_mpoly_mul_monomial(k2, f, g, mctx, cctx);

            if (status == GR_SUCCESS)
                gr_mpoly_assert_canonical(k2, mctx, cctx);

            status |= gr_mpoly_mul_monomial(t, f, h, mctx, cctx);

            if (status == GR_SUCCESS)
                gr_mpoly_assert_canonical(t, mctx, cctx);

            status |= gr_mpoly_add(k2, k2, t, mctx, cctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(k1, k2, mctx, cctx) == T_FALSE)
            {
                flint_printf("FAIL: Check (f + g) - g = f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("h = "); gr_mpoly_print_pretty(h, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f * (g + h) = "); gr_mpoly_print_pretty(k1, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f * g + f * h = "); gr_mpoly_print_pretty(k2, NULL, mctx, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_mpoly_clear(g, mctx, cctx);
        gr_mpoly_clear(h, mctx, cctx);
        gr_mpoly_clear(k1, mctx, cctx);
        gr_mpoly_clear(k2, mctx, cctx);
        gr_mpoly_clear(t, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
