/*
    Copyright (C) 2026

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_mpoly_sqr, state)
{
    slong i, j;

    gr_ctx_t ctx1;
    gr_ctx_init_nmod(ctx1, 17);

    /* Check gr_mpoly_sqr(f) == gr_mpoly_mul_johnson(f, f) for every ring:
       for (approx-)commutative rings the dispatcher routes to the dedicated
       squaring code, for noncommutative rings it falls through to the general
       product; either way the identity must hold. */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, h1, h2;
        slong len1;
        flint_bitcnt_t exp_bits1;
        int status, aliasing;

        switch (n_randint(state, 6))
        {
            case 0:
                gr_ctx_init_random(cctx, state);            /* anything */
                break;
            case 1:
                gr_ctx_init_fmpz(cctx);                     /* exact commutative */
                break;
            case 2:
                gr_ctx_init_real_arb(cctx, 2 + n_randint(state, 200)); /* approx */
                break;
            case 3:
                gr_ctx_init_complex_acb(cctx, 2 + n_randint(state, 200)); /* approx */
                break;
            case 4:
                gr_ctx_init_random_finite_field(cctx, state);
                break;
            default:
                gr_ctx_init_matrix_ring(cctx, ctx1, 2);     /* noncommutative */
                break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(h1, ctx);
        gr_mpoly_init(h2, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        /* mix of lengths so we hit the monomial, johnson and threaded paths */
        if (gr_ctx_is_finite(cctx) == T_TRUE)
        {
            switch (n_randint(state, 4))
            {
                case 0:  len1 = n_randint(state, 2); break;   /* 0 or 1 -> monomial */
                case 1:  len1 = n_randint(state, 20); break;
                default: len1 = n_randint(state, 400); break;
            }
        }
        else
        {
            len1 = n_randint(state, 10);
        }

        exp_bits1 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(h1, state, n_randint(state, 10),
                                                   n_randint(state, 100) + 2, ctx);
            status |= gr_mpoly_randtest_bits(h2, state, n_randint(state, 10),
                                                   n_randint(state, 100) + 2, ctx);

            /* reference */
            status |= gr_mpoly_mul_johnson(h1, f, f, ctx);

            aliasing = n_randint(state, 2);

            if (aliasing == 0)
            {
                status |= gr_mpoly_sqr(h2, f, ctx);
            }
            else
            {
                status |= gr_mpoly_set(h2, f, ctx);
                status |= gr_mpoly_sqr(h2, h2, ctx);        /* in place */
            }

            if (status == GR_SUCCESS)
            {
                gr_mpoly_assert_canonical(h1, ctx);
                gr_mpoly_assert_canonical(h2, ctx);

                if (gr_mpoly_equal(h1, h2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: sqr(f) != f * f\n");
                    flint_printf("i = %wd, j = %wd, aliasing = %d\n", i, j, aliasing);
                    gr_ctx_println(ctx);
                    flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                    flint_printf("f*f    = "); gr_mpoly_print_pretty(h1, ctx); flint_printf("\n");
                    flint_printf("sqr(f) = "); gr_mpoly_print_pretty(h2, ctx); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(h1, ctx);
        gr_mpoly_clear(h2, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Focused check of gr_mpoly_sqr_monomial, including in-place aliasing,
       against the general product. */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, h1, h2;
        int status;

        gr_ctx_init_random(cctx, state);
        gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(h1, ctx);
        gr_mpoly_init(h2, ctx);

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            /* a single-term polynomial (or empty) */
            status |= gr_mpoly_randtest_bits(f, state, 1,
                                              n_randint(state, 100) + 2, ctx);

            status |= gr_mpoly_mul_johnson(h1, f, f, ctx);

            /* out-of-place */
            status |= gr_mpoly_sqr_monomial(h2, f, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(h1, h2, ctx) == T_FALSE)
            {
                flint_printf("FAIL: sqr_monomial out-of-place\n");
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            /* in-place */
            status |= gr_mpoly_set(h2, f, ctx);
            status |= gr_mpoly_sqr_monomial(h2, h2, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(h1, h2, ctx) == T_FALSE)
            {
                flint_printf("FAIL: sqr_monomial in-place\n");
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(h1, ctx);
        gr_mpoly_clear(h2, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    gr_ctx_clear(ctx1);

    TEST_FUNCTION_END(state);
}
