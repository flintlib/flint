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

TEST_FUNCTION_START(gr_mpoly_divides_heap_threaded, state)
{
    slong i, j;

    gr_ctx_t ctx1;
    gr_ctx_init_nmod(ctx1, 17);

    /* Check that the threaded quotient matches the serial heap quotient */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_ctx_t debug_base_ctx;
        int used_debug_base = 0;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, h, q1, q2;
        slong lenf, leng;
        flint_bitcnt_t exp_bits;
        int status, aliasing;
        int big = (n_randint(state, 8) == 0);

        switch (n_randint(state, 5))
        {
            case 0:
                gr_ctx_init_random_finite_field(cctx, state);
                break;
            case 1:
                gr_ctx_init_fmpz(cctx);
                break;
            case 2:
                gr_ctx_init_fmpq(cctx);
                break;
            case 3:
                /* Regression coverage: a ring whose operations can
                   randomly return GR_UNABLE must never cause a spurious
                   GR_DOMAIN verdict on an exactly-divisible instance (the
                   "not exact" detection must OR GR_DOMAIN into any status
                   already accumulated while computing the tainted term,
                   not overwrite it -- see gr_mpoly_divides_heap). Keep
                   sizes modest here: GR_DEBUG_WRAP gives every coefficient
                   a real heap allocation, making "big" instances far
                   heavier per-term than for the other rings tested above. */
                gr_ctx_init_nmod8(debug_base_ctx, 251);
                used_debug_base = 1;
                gr_ctx_init_debug(cctx, debug_base_ctx, n_randint(state, 2),
                                   n_randint(state, 2) ? 0.0 : 0.001);
                big = 0;
                break;
            default:
                gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1));
                break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(h, ctx);
        gr_mpoly_init(q1, ctx);
        gr_mpoly_init(q2, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        lenf = n_randint(state, big ? 80 : 20) + 1;
        leng = n_randint(state, big ? 30 : 10) + 1;
        exp_bits = n_randint(state, big ? 12 : 4) + 2;

        for (j = 0; j < 2; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, lenf, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(g, state, leng, exp_bits, ctx);

            if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
                status |= gr_mpoly_one(g, ctx);

            /* h = f * g, so g is guaranteed to divide h with quotient f */
            status |= gr_mpoly_mul(h, f, g, ctx);

            if (status != GR_SUCCESS)
                continue;

            status = GR_SUCCESS;
            status |= gr_mpoly_randtest_bits(q1, state, n_randint(state, 3),
                                                   n_randint(state, 4) + 2, ctx);
            status |= gr_mpoly_randtest_bits(q2, state, n_randint(state, 3),
                                                   n_randint(state, 4) + 2, ctx);

            {
                int dstatus1 = gr_mpoly_divides_heap(q1, h, g, ctx);
                int dstatus2;

                aliasing = n_randint(state, 3);

                switch (aliasing)
                {
                    case 0:
                        dstatus2 = gr_mpoly_divides_heap_threaded(q2, h, g, ctx);
                        break;
                    case 1:
                        status |= gr_mpoly_set(q2, h, ctx);
                        dstatus2 = gr_mpoly_divides_heap_threaded(q2, q2, g, ctx);
                        break;
                    default:
                        status |= gr_mpoly_set(q2, g, ctx);
                        dstatus2 = gr_mpoly_divides_heap_threaded(q2, h, q2, ctx);
                        break;
                }

                if (status != GR_SUCCESS)
                    continue;

                if (dstatus1 == GR_SUCCESS && dstatus2 == GR_SUCCESS)
                {
                    gr_mpoly_assert_canonical(q1, ctx);
                    gr_mpoly_assert_canonical(q2, ctx);

                    if (gr_mpoly_equal(q1, q2, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: threaded != heap\n");
                        flint_printf("i = %wd, j = %wd, aliasing = %d\n", i, j, aliasing);
                        gr_ctx_println(cctx);
                        flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                        flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                        flint_printf("heap     = "); gr_mpoly_print_pretty(q1, ctx); flint_printf("\n");
                        flint_printf("threaded = "); gr_mpoly_print_pretty(q2, ctx); flint_printf("\n");
                        fflush(stdout);
                        flint_abort();
                    }

                    if (gr_mpoly_equal(q1, f, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: quotient != f\n");
                        flint_printf("i = %wd, j = %wd\n", i, j);
                        gr_ctx_println(cctx);
                        fflush(stdout);
                        flint_abort();
                    }
                }
                else if (dstatus1 == GR_DOMAIN || dstatus2 == GR_DOMAIN)
                {
                    /* g always divides h = f*g exactly; a definite GR_DOMAIN
                       verdict from either algorithm would be a bug */
                    flint_printf("FAIL: exact division incorrectly reported as GR_DOMAIN\n");
                    flint_printf("i = %wd, j = %wd, dstatus1 = %d, dstatus2 = %d\n",
                                                            i, j, dstatus1, dstatus2);
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(h, ctx);
        gr_mpoly_clear(q1, ctx);
        gr_mpoly_clear(q2, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
        if (used_debug_base)
            gr_ctx_clear(debug_base_ctx);
    }

    /* Also check on arbitrary (not necessarily divisible) inputs that the
       two algorithms never disagree when both succeed */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, q1, q2;
        int status, dstatus1, dstatus2;
        int big = (n_randint(state, 8) == 0);

        gr_mpoly_ctx_init_rand(ctx, state, ctx1, 4);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(q1, ctx);
        gr_mpoly_init(q2, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(f, state, n_randint(state, big ? 100 : 20) + 1,
                                                n_randint(state, big ? 12 : 4) + 2, ctx);
        status |= gr_mpoly_randtest_bits(g, state, n_randint(state, big ? 30 : 10) + 1,
                                                n_randint(state, big ? 12 : 4) + 2, ctx);
        if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
            status |= gr_mpoly_one(g, ctx);

        if (status == GR_SUCCESS)
        {
            dstatus1 = gr_mpoly_divides_heap(q1, f, g, ctx);
            dstatus2 = gr_mpoly_divides_heap_threaded(q2, f, g, ctx);

            if (dstatus1 != dstatus2)
            {
                flint_printf("FAIL: status mismatch (heap=%d, threaded=%d)\n",
                                                            dstatus1, dstatus2);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (dstatus1 == GR_SUCCESS)
            {
                gr_mpoly_assert_canonical(q1, ctx);
                gr_mpoly_assert_canonical(q2, ctx);

                if (gr_mpoly_equal(q1, q2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: threaded != heap (arbitrary inputs)\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(q1, ctx);
        gr_mpoly_clear(q2, ctx);

        gr_mpoly_ctx_clear(ctx);
    }

    gr_ctx_clear(ctx1);

    TEST_FUNCTION_END(state);
}
