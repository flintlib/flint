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

TEST_FUNCTION_START(gr_mpoly_divrem, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t a, b, q, r, t, u;
        slong len, blen;
        flint_bitcnt_t ebits;
        int status, big;
        truth_t is_field;

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
        else
            gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1));

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
                status = gr_mpoly_divrem(q, r, a, b, ctx);
                break;
            case 1:
                status = gr_mpoly_set(q, b, ctx);
                status |= gr_mpoly_divrem(q, r, a, q, ctx);
                break;
            case 2:
                status = gr_mpoly_set(r, a, ctx);
                status |= gr_mpoly_divrem(q, r, r, b, ctx);
                break;
            default:
                status = gr_mpoly_set(r, a, ctx);
                status = gr_mpoly_set(q, b, ctx);
                status |= gr_mpoly_divrem(q, r, r, q, ctx);
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
            if (gr_mpoly_div(t, a, b, ctx) == GR_SUCCESS &&
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

        /* nonfield (euclidean) variant: a == q*b + r always */
        status = gr_mpoly_divrem_weak(q, r, a, b, ctx);
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

            if (gr_mpoly_div_weak(t, a, b, ctx) == GR_SUCCESS &&
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
    }

    TEST_FUNCTION_END(state);
}
