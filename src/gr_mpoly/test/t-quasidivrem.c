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

TEST_FUNCTION_START(gr_mpoly_quasidivrem, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t a, b, q, r, t, u;
        gr_ptr scale, scale2;
        slong len;
        flint_bitcnt_t ebits;
        int status;

        /* fmpz exercises the scaling; the others are fields (scale == 1) */
        switch (n_randint(state, 4))
        {
            case 0: gr_ctx_init_fmpz(cctx); break;
            case 1: gr_ctx_init_fmpq(cctx); break;
            case 2: gr_ctx_init_random_finite_field(cctx, state); break;
            default: gr_ctx_init_nmod(cctx, n_randtest_prime(state, 1)); break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 12);

        gr_mpoly_init(a, ctx);
        gr_mpoly_init(b, ctx);
        gr_mpoly_init(q, ctx);
        gr_mpoly_init(r, ctx);
        gr_mpoly_init(t, ctx);
        gr_mpoly_init(u, ctx);
        GR_TMP_INIT(scale, cctx);
        GR_TMP_INIT(scale2, cctx);

        /* keep inputs small: over Z the scale grows like lc(B)^(#steps) */
        len = 1 + n_randint(state, 7);
        ebits = 2 + n_randint(state, 5);

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(a, state, len, ebits, ctx);
        status |= gr_mpoly_randtest_bits(b, state, 1 + n_randint(state, 6), ebits, ctx);
        if (gr_mpoly_is_zero(b, ctx) != T_FALSE)
            status |= gr_mpoly_one(b, ctx);

        if (status != GR_SUCCESS)
            goto next;

        /* field-like */
        switch (n_randint(state, 4))
        {
            case 0:
                status = gr_mpoly_quasidivrem(scale, q, r, a, b, ctx);
                break;
            case 1:
                status = gr_mpoly_set(q, b, ctx);
                status |= gr_mpoly_quasidivrem(scale, q, r, a, q, ctx);
                break;
            case 2:
                status = gr_mpoly_set(r, a, ctx);
                status |= gr_mpoly_quasidivrem(scale, q, r, r, b, ctx);
                break;
            default:
                status = gr_mpoly_set(r, a, ctx);
                status = gr_mpoly_set(q, b, ctx);
                status |= gr_mpoly_quasidivrem(scale, q, r, r, q, ctx);
                break;
        }

        if (status == GR_SUCCESS)
        {
            gr_mpoly_assert_canonical(q, ctx);
            gr_mpoly_assert_canonical(r, ctx);

            /* check scale*a == q*b + r */
            if (gr_mpoly_mul(t, q, b, ctx) == GR_SUCCESS &&
                gr_mpoly_add(t, t, r, ctx) == GR_SUCCESS &&
                gr_mpoly_mul_scalar(u, a, scale, ctx) == GR_SUCCESS &&
                gr_mpoly_equal(t, u, ctx) == T_FALSE)
            {
                flint_printf("FAIL: scale*a != q*b + r\n");
                gr_ctx_println(cctx);
                flint_printf("a = "); gr_mpoly_print_pretty(a, ctx); flint_printf("\n");
                flint_printf("b = "); gr_mpoly_print_pretty(b, ctx); flint_printf("\n");
                flint_printf("q = "); gr_mpoly_print_pretty(q, ctx); flint_printf("\n");
                flint_printf("r = "); gr_mpoly_print_pretty(r, ctx); flint_printf("\n");
                flint_printf("scale = "); gr_println(scale, cctx);
                fflush(stdout);
                flint_abort();
            }

            /* quasidiv must agree on the quotient and the scale */
            if (gr_mpoly_quasidiv(scale2, t, a, b, ctx) == GR_SUCCESS)
            {
                if (gr_mpoly_equal(t, q, ctx) == T_FALSE ||
                    gr_equal(scale2, scale, cctx) == T_FALSE)
                {
                    flint_printf("FAIL: quasidiv disagrees with quasidivrem\n");
                    gr_ctx_println(cctx);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

next:
        gr_mpoly_clear(a, ctx);
        gr_mpoly_clear(b, ctx);
        gr_mpoly_clear(q, ctx);
        gr_mpoly_clear(r, ctx);
        gr_mpoly_clear(t, ctx);
        gr_mpoly_clear(u, ctx);
        GR_TMP_CLEAR(scale, cctx);
        GR_TMP_CLEAR(scale2, cctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
