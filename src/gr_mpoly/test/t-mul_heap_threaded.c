/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_mpoly_mul_heap_threaded, state)
{
    slong i, j;

    gr_ctx_t ctx1;
    gr_ctx_init_nmod(ctx1, 17);

    /* Check that the threaded product matches the serial Johnson product */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, h1, h2;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;
        int status, aliasing;

        switch (n_randint(state, 4))
        {
            case 0:
                gr_ctx_init_random_finite_field(cctx, state);
                break;
            case 1:
                gr_ctx_init_fmpz(cctx);
                break;
            case 2:
                gr_ctx_init_real_arb(cctx, 2 + n_randint(state, 200));
                break;
            case 3:   /* noncommutative ring */
                gr_ctx_init_matrix_ring(cctx, ctx1, 2);
                break;

        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(h1, ctx);
        gr_mpoly_init(h2, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        len1 = n_randint(state, 500);
        len2 = n_randint(state, 200);

        if (n_randint(state, 4) == 0)
        {
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 10) + 2;
        }
        else
        {
            exp_bits1 = n_randint(state, 200) + 2;
            exp_bits2 = n_randint(state, 200) + 2;
        }

        for (j = 0; j < 4; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h1, state, n_randint(state, 10),
                                                   n_randint(state, 100) + 2, ctx);
            status |= gr_mpoly_randtest_bits(h2, state, n_randint(state, 10),
                                                   n_randint(state, 100) + 2, ctx);

            status |= gr_mpoly_mul_heap(h1, f, g, ctx);

            aliasing = n_randint(state, 3);

            switch (aliasing)
            {
                case 0:
                    status |= gr_mpoly_mul_heap_threaded(h2, f, g, ctx);
                    break;
                case 1:
                    status |= gr_mpoly_set(h2, f, ctx);
                    status |= gr_mpoly_mul_heap_threaded(h2, h2, g, ctx);
                    break;
                default:
                    status |= gr_mpoly_set(h2, g, ctx);
                    status |= gr_mpoly_mul_heap_threaded(h2, f, h2, ctx);
                    break;
            }

            if (status == GR_SUCCESS)
            {
                gr_mpoly_assert_canonical(h1, ctx);
                gr_mpoly_assert_canonical(h2, ctx);

                if (gr_mpoly_equal(h1, h2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: threaded != johnson\n");
                    flint_printf("i = %wd, j = %wd, aliasing = %d\n", i, j, aliasing);
                    gr_ctx_println(ctx);
                    flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                    flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                    flint_printf("johnson  = "); gr_mpoly_print_pretty(h1, ctx); flint_printf("\n");
                    flint_printf("threaded = "); gr_mpoly_print_pretty(h2, ctx); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(h1, ctx);
        gr_mpoly_clear(h2, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    gr_ctx_clear(ctx1);

    TEST_FUNCTION_END(state);
}
