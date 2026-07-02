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

static void
_check_one(flint_rand_t state, gr_ctx_t cctx, gr_mpoly_ctx_t ctx, int big)
{
    gr_mpoly_t f, g, h, q, t;
    slong len, lenf, leng;
    flint_bitcnt_t ebits;
    int status, dstatus;
    truth_t dom;

    gr_mpoly_init(f, ctx);
    gr_mpoly_init(g, ctx);
    gr_mpoly_init(h, ctx);
    gr_mpoly_init(q, ctx);
    gr_mpoly_init(t, ctx);

    dom = gr_mpoly_ctx_is_integral_domain(ctx);

    len = big ? 30 : 8;
    lenf = n_randint(state, len) + 1;
    leng = n_randint(state, len) + 1;
    ebits = 2 + n_randint(state, big ? 10 : 3);

    status = GR_SUCCESS;
    status |= gr_mpoly_randtest_bits(f, state, lenf, ebits, ctx);
    status |= gr_mpoly_randtest_bits(g, state, leng, ebits, ctx);

    /* make the divisor nonzero */
    if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
        status |= gr_mpoly_one(g, ctx);

    /* h = f * g is divisible by g with quotient f */
    status |= gr_mpoly_mul(h, f, g, ctx);

    if (status == GR_SUCCESS)
    {
        dstatus = gr_mpoly_divides(q, h, g, ctx);

        /* over an integral domain, g | f*g must hold (unless arithmetic
           was unable to decide) */
        if (dom == T_TRUE && dstatus == GR_DOMAIN)
        {
            flint_printf("FAIL: integral domain, g should divide f*g\n");
            gr_ctx_println(cctx);
            flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
            flint_printf("h = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        /* whenever divisibility is reported, q*g must reproduce h */
        if (dstatus == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q, ctx);

            if (gr_mpoly_mul(t, q, g, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, h, ctx) == T_FALSE)
            {
                flint_printf("FAIL: q*g != h\n");
                gr_ctx_println(cctx);
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("h = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
                flint_printf("q = "); gr_mpoly_print_pretty(q, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* test aliasing of quotient with dividend */
        dstatus = gr_mpoly_divides(h, h, g, ctx);
        if (dstatus == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(h, ctx);
            if (gr_mpoly_equal(h, q, ctx) == T_FALSE)
            {
                flint_printf("FAIL: aliasing q == h/g\n");
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    /* consistency on arbitrary inputs: if divisible, q*g == f */
    status = GR_SUCCESS;
    status |= gr_mpoly_randtest_bits(f, state, n_randint(state, len) + 1, ebits, ctx);
    status |= gr_mpoly_randtest_bits(g, state, n_randint(state, len) + 1, ebits, ctx);
    if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
        status |= gr_mpoly_one(g, ctx);

    if (status == GR_SUCCESS)
    {
        dstatus = gr_mpoly_divides(q, f, g, ctx);
        if (dstatus == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q, ctx);
            if (gr_mpoly_mul(t, q, g, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, f, ctx) == T_FALSE)
            {
                flint_printf("FAIL: reported divisible but q*g != f\n");
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("q = "); gr_mpoly_print_pretty(q, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    /* explicitly exercise the scalar/monomial divisor special case */
    status = GR_SUCCESS;

    if (n_randint(state, 4) == 0)
        ebits = 2 + n_randint(state, 200);
    status |= gr_mpoly_randtest_bits(f, state, n_randint(state, len) + 1, ebits, ctx);
    status |= gr_mpoly_randtest_bits(g, state, 1, ebits, ctx);   /* monomial or scalar */
    if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
        status |= gr_mpoly_one(g, ctx);
    status |= gr_mpoly_mul(h, f, g, ctx);

    if (status == GR_SUCCESS)
    {
        dstatus = gr_mpoly_divides(q, h, g, ctx);

        if (dom == T_TRUE && dstatus == GR_DOMAIN)
        {
            flint_printf("FAIL: monomial divisor should divide f*g\n");
            gr_ctx_println(cctx);
            fflush(stdout);
            flint_abort();
        }

        if (dstatus == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q, ctx);
            if (gr_mpoly_mul(t, q, g, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, h, ctx) == T_FALSE)
            {
                flint_printf("FAIL: monomial divisor q*g != h\n");
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }
        }

        /* arbitrary dividend by the same monomial */
        dstatus = gr_mpoly_divides(q, f, g, ctx);
        if (dstatus == T_TRUE)
        {
            gr_mpoly_assert_canonical(q, ctx);
            if (gr_mpoly_mul(t, q, g, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, f, ctx) == T_FALSE)
            {
                flint_printf("FAIL: monomial divisor reported divisible but q*g != f\n");
                gr_ctx_println(cctx);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    gr_mpoly_clear(f, ctx);
    gr_mpoly_clear(g, ctx);
    gr_mpoly_clear(h, ctx);
    gr_mpoly_clear(q, ctx);
    gr_mpoly_clear(t, ctx);
}

TEST_FUNCTION_START(gr_mpoly_divides, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        int big = (n_randint(state, 8) == 0);

        /* commutative rings only: choose a field or integral domain */
        switch (n_randint(state, 5))
        {
            case 0:
                gr_ctx_init_fmpz(cctx);
                break;
            case 1:
                gr_ctx_init_fmpq(cctx);
                break;
            case 2:
                gr_ctx_init_random_finite_field(cctx, state);
                break;
            case 3:
                gr_ctx_init_real_arb(cctx, 2 + n_randint(state, 200));
                break;
            default:
                gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1));
                break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 20);

        flint_set_num_threads(1 + n_randint(state, 4));

        _check_one(state, cctx, ctx, big);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
