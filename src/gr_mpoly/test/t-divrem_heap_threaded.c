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
    slong i, iter;

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

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t a, b, q, r, t, u;
        slong len, blen;
        flint_bitcnt_t ebits;
        int status, big;
        truth_t is_field;

        flint_set_num_threads(1 + n_randint(state, 4));

        gr_ctx_t ctx1;
        gr_ctx_init_nmod8(ctx1, n_randprime(state, 8, 1));

        big = (n_randint(state, 4) == 0);

        if (big)
        {
            /* bounded-coefficient ring + large exponents (few terms) to
               exercise the overflow retry without coefficient blow-up */
            gr_ctx_init_nmod(cctx, n_randprime(state, 5 + n_randint(state, 10), 1));
        }
        else if (n_randint(state, 2))
            gr_ctx_init_random_finite_field(cctx, state);
        else if (n_randint(state, 3) == 0)
            gr_ctx_init_fmpz(cctx);
        else if (n_randint(state, 2))
            gr_ctx_init_fmpq(cctx);
        else if (n_randint(state, 2))
            gr_ctx_init_real_arb(cctx, 2 + n_randint(state, 200));
        else if (n_randint(state, 2))
            gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1));
        else
            gr_ctx_init_debug(cctx, ctx1, n_randint(state, 2), n_randint(state, 2) ? 0.0 : 0.001);

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 15);
        is_field = gr_mpoly_ctx_is_field(ctx);

        gr_mpoly_init(a, ctx);
        gr_mpoly_init(b, ctx);
        gr_mpoly_init(q, ctx);
        gr_mpoly_init(r, ctx);
        gr_mpoly_init(t, ctx);
        gr_mpoly_init(u, ctx);

        if (big)
        {
            len = 1 + n_randint(state, 4);
            blen = 1 + n_randint(state, 4);
            ebits = 2 + n_randint(state, 10);
        }
        else
        {
            len = 1 + n_randint(state, 12);
            blen = 1 + n_randint(state, 12);
            ebits = 2 + n_randint(state, 4);
        }

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(a, state, len, ebits, ctx);
        status |= gr_mpoly_randtest_bits(b, state, blen, ebits, ctx);

        if (gr_mpoly_is_zero(b, ctx) != T_FALSE)
            status |= gr_mpoly_one(b, ctx);

        if (status != GR_SUCCESS)
            goto next;

        switch (n_randint(state, 4))
        {
            case 0:
                status = gr_mpoly_divrem_heap_threaded(q, r, a, b, ctx);
                break;
            case 1:
                status = gr_mpoly_set(q, b, ctx);
                status |= gr_mpoly_divrem_heap_threaded(q, r, a, q, ctx);
                break;
            case 2:
                status = gr_mpoly_set(r, a, ctx);
                status |= gr_mpoly_divrem_heap_threaded(q, r, r, b, ctx);
                break;
            default:
                status = gr_mpoly_set(r, a, ctx);
                status = gr_mpoly_set(q, b, ctx);
                status |= gr_mpoly_divrem_heap_threaded(q, r, r, q, ctx);
                break;
        }

        if (status == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q, ctx);
            gr_mpoly_assert_canonical(r, ctx);

            /* check a == q*b + r */
            if (gr_mpoly_mul(t, q, b, ctx) == GR_SUCCESS &&
                gr_mpoly_add(t, t, r, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, a, ctx) == T_FALSE)
            {
                flint_printf("FAIL: a != q*b + r\n");
                gr_ctx_println(cctx);
                flint_printf("a = "); gr_mpoly_print_pretty(a, ctx); flint_printf("\n");
                flint_printf("b = "); gr_mpoly_print_pretty(b, ctx); flint_printf("\n");
                flint_printf("q = "); gr_mpoly_print_pretty(q, ctx); flint_printf("\n");
                flint_printf("r = "); gr_mpoly_print_pretty(r, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            /* div must agree with divrem on the quotient */
            if (gr_mpoly_div_heap_threaded(t, a, b, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, q, ctx) == T_FALSE)
            {
                flint_printf("FAIL: div != divrem quotient\n");
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }
        }
        else if (is_field == T_TRUE && status == GR_DOMAIN)
        {
            flint_printf("FAIL: divrem returned GR_DOMAIN over a field\n");
            gr_ctx_println(cctx);
            fflush(stdout);
            flint_abort();
        }

        if (cctx->which_ring != GR_CTX_DEBUG &&
            status != GR_SUCCESS &&
            gr_is_invertible(b->coeffs, cctx) == T_TRUE)
        {
            flint_printf("FAIL: division unable although lc(g) is invertible\n");
            gr_ctx_println(ctx);
            fflush(stdout);
            flint_abort();
        }

        /* nonfield (euclidean) variant: a == q*b + r always */
        status = gr_mpoly_divrem_weak_heap_threaded(q, r, a, b, ctx);
        if (status == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q, ctx);
            gr_mpoly_assert_canonical(r, ctx);

            if (gr_mpoly_mul(t, q, b, ctx) == GR_SUCCESS &&
                gr_mpoly_add(t, t, r, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, a, ctx) == T_FALSE)
            {
                flint_printf("FAIL: weak a != q*b + r\n");
                gr_ctx_println(cctx);
                flint_printf("a = "); gr_mpoly_print_pretty(a, ctx); flint_printf("\n");
                flint_printf("b = "); gr_mpoly_print_pretty(b, ctx); flint_printf("\n");
                flint_printf("q = "); gr_mpoly_print_pretty(q, ctx); flint_printf("\n");
                flint_printf("r = "); gr_mpoly_print_pretty(r, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (gr_mpoly_div_weak_heap_threaded(t, a, b, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, q, ctx) == T_FALSE)
            {
                flint_printf("FAIL: div_weak != divrem_weak quotient\n");
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }

            /* over a field the two variants coincide */
            if (is_field == T_TRUE &&
                gr_mpoly_divrem(t, u, a, b, ctx) == GR_SUCCESS &&
                (gr_mpoly_equal(t, q, ctx) == T_FALSE ||
                 gr_mpoly_equal(u, r, ctx) == T_FALSE))
            {
                flint_printf("FAIL: field/weak disagree over a field\n");
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }
        }

next:
        gr_mpoly_clear(a, ctx);
        gr_mpoly_clear(b, ctx);
        gr_mpoly_clear(q, ctx);
        gr_mpoly_clear(r, ctx);
        gr_mpoly_clear(t, ctx);
        gr_mpoly_clear(u, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
        gr_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
