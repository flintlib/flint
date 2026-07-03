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

TEST_FUNCTION_START(gr_mpoly_divrem_heap_threaded, state)
{
    slong i;

    /* Check that the threaded div/divrem family matches the serial kernel,
       for both the exact-coefficient and Euclidean ("weak") variants, with
       and without a remainder. */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_ctx_t debug_base_ctx;
        int used_debug_base = 0;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, q1, q2, r1, r2, t;
        slong lenf, leng;
        flint_bitcnt_t exp_bits;
        int status, nonfield, want_r;
        int s1, s2;
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
                   GR_DOMAIN verdict on an exactly-divisible instance --
                   see the analogous case in t-divides_heap_threaded.c and
                   the "status |= cstatus" pattern in
                   STRIPE_DIVREM_COEFF_STEP / DIVREM_COEFF_STEP. Keep sizes
                   modest: GR_DEBUG_WRAP gives every coefficient a real
                   heap allocation, making "big" instances far heavier per
                   term than for the other rings tested above. */
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
        gr_mpoly_init(q1, ctx);
        gr_mpoly_init(q2, ctx);
        gr_mpoly_init(r1, ctx);
        gr_mpoly_init(r2, ctx);
        gr_mpoly_init(t, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        lenf = n_randint(state, big ? 80 : 20) + 1;
        leng = n_randint(state, big ? 30 : 10) + 1;
        exp_bits = n_randint(state, big ? 12 : 4) + 2;

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(f, state, lenf, exp_bits, ctx);
        status |= gr_mpoly_randtest_bits(g, state, leng, exp_bits, ctx);
        if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
            status |= gr_mpoly_one(g, ctx);

        if (status != GR_SUCCESS)
        {
            gr_mpoly_clear(f, ctx); gr_mpoly_clear(g, ctx);
            gr_mpoly_clear(q1, ctx); gr_mpoly_clear(q2, ctx);
            gr_mpoly_clear(r1, ctx); gr_mpoly_clear(r2, ctx);
            gr_mpoly_clear(t, ctx);
            gr_mpoly_ctx_clear(ctx); gr_ctx_clear(cctx);
            if (used_debug_base) gr_ctx_clear(debug_base_ctx);
            continue;
        }

        nonfield = n_randint(state, 2);
        want_r = n_randint(state, 2);

        if (want_r)
        {
            s1 = nonfield ? gr_mpoly_divrem_weak(q1, r1, f, g, ctx)
                          : gr_mpoly_divrem_heap(q1, r1, f, g, ctx);
            s2 = nonfield ? gr_mpoly_divrem_weak_heap_threaded(q2, r2, f, g, ctx)
                          : gr_mpoly_divrem_heap_threaded(q2, r2, f, g, ctx);
        }
        else
        {
            s1 = nonfield ? gr_mpoly_div_weak(q1, f, g, ctx)
                          : gr_mpoly_div(q1, f, g, ctx);
            s2 = nonfield ? gr_mpoly_div_weak_heap_threaded(q2, f, g, ctx)
                          : gr_mpoly_div_heap_threaded(q2, f, g, ctx);
        }

        if (s1 != s2)
        {
            flint_printf("FAIL: status mismatch serial=%d threaded=%d "
                         "nonfield=%d want_r=%d\n", s1, s2, nonfield, want_r);
            flint_printf("i = %wd\n", i);
            gr_ctx_println(cctx);
            flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if (s1 == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q1, ctx);
            gr_mpoly_assert_canonical(q2, ctx);

            if (gr_mpoly_equal(q1, q2, ctx) == T_FALSE)
            {
                flint_printf("FAIL: quotient mismatch nonfield=%d want_r=%d\n",
                                                            nonfield, want_r);
                flint_printf("i = %wd\n", i);
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }

            if (want_r)
            {
                gr_mpoly_assert_canonical(r1, ctx);
                gr_mpoly_assert_canonical(r2, ctx);

                if (gr_mpoly_equal(r1, r2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: remainder mismatch nonfield=%d\n", nonfield);
                    flint_printf("i = %wd\n", i);
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }

                /* reconstruction: q2*g + r2 == f */
                if (gr_mpoly_mul(t, q2, g, ctx) == GR_SUCCESS &&
                    gr_mpoly_add(t, t, r2, ctx) == GR_SUCCESS)
                {
                    if (gr_mpoly_equal(t, f, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: q*g + r != f\n");
                        flint_printf("i = %wd\n", i);
                        gr_ctx_println(cctx);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(q1, ctx);
        gr_mpoly_clear(q2, ctx);
        gr_mpoly_clear(r1, ctx);
        gr_mpoly_clear(r2, ctx);
        gr_mpoly_clear(t, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
        if (used_debug_base)
            gr_ctx_clear(debug_base_ctx);
    }

    TEST_FUNCTION_END(state);
}
