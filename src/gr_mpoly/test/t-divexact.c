/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_divexact, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, h, q, r;
        slong lenf, leng;
        flint_bitcnt_t exp_bits;
        int status, s1;
        int monomial_case = (n_randint(state, 4) == 0);
        int big = (!monomial_case && n_randint(state, 8) == 0);

        gr_ctx_t ctx1;
        gr_ctx_init_nmod8(ctx1, n_randprime(state, 8, 1));

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
                gr_ctx_init_real_arb(cctx, 2 + n_randint(state, 200));
                break;
            default:
                gr_ctx_init_debug(cctx, ctx1, n_randint(state, 2), n_randint(state, 2) ? 0.0 : 0.001);
                break;
        }

        gr_mpoly_ctx_init_rand(ctx, state, cctx, 12);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(h, ctx);
        gr_mpoly_init(q, ctx);
        gr_mpoly_init(r, ctx);

        flint_set_num_threads(1 + n_randint(state, 4));

        /* occasionally big enough to actually exercise the threaded
           dispatch inside gr_mpoly_divexact (A->length > 500, B->length
           > 2), and always call gr_mpoly_divexact_heap_threaded directly
           too, regardless of size, to exercise it even when small inputs
           would otherwise make gr_mpoly_divexact stay serial */
        lenf = n_randint(state, big ? 400 : 20) + 1;
        leng = monomial_case ? 1 : (n_randint(state, big ? 80 : 10) + 3);
        exp_bits = n_randint(state, big ? 8 : 4) + 2;

        status = GR_SUCCESS;
        status |= gr_mpoly_randtest_bits(f, state, lenf, exp_bits, ctx);
        status |= gr_mpoly_randtest_bits(g, state, leng, exp_bits, ctx);
        if (gr_mpoly_is_zero(g, ctx) != T_FALSE)
            status |= gr_mpoly_one(g, ctx);

        if (status == GR_SUCCESS)
        {
            /* h = f * g, so g divides h exactly with quotient f */
            status = gr_mpoly_mul(h, f, g, ctx);

            if (status == GR_SUCCESS)
            {
                int alg;

                for (alg = 0; alg < 6; alg++)
                {
                    if (alg == 0)
                        s1 = gr_mpoly_divides(q, h, g, ctx);
                    else if (alg == 1)
                        s1 = gr_mpoly_div(q, h, g, ctx);
                    else if (alg == 2)
                        s1 = gr_mpoly_divexact(q, h, g, ctx);
                    else if (alg == 3)
                        s1 = gr_mpoly_divexact_heap_threaded(q, h, g, ctx);
                    else if (alg == 4)
                        s1 = gr_mpoly_divrem(q, r, h, g, ctx);
                    else if (alg == 5)
                        s1 = gr_mpoly_divrem_heap_threaded(q, r, h, g, ctx);

                    gr_mpoly_assert_canonical(q, ctx);

                    if (s1 == GR_DOMAIN)
                    {
                        /* g always divides h = f*g exactly */
                        flint_printf("FAIL: exact division reported GR_DOMAIN\n");
                        flint_printf("alg = %d\n", alg);
                        gr_ctx_println(ctx);
                        fflush(stdout);
                        flint_abort();
                    }

                    if (cctx->which_ring != GR_CTX_DEBUG &&
                        !(alg == 0 && cctx->which_ring == GR_CTX_RR_ARB) &&
                        s1 != GR_SUCCESS &&
                        gr_is_invertible(g->coeffs, cctx) == T_TRUE)
                    {
                        flint_printf("FAIL: division unable although lc(g) is invertible\n");
                        flint_printf("alg = %d\n", alg);
                        gr_ctx_println(ctx);
                        fflush(stdout);
                        flint_abort();
                    }

                    if (s1 == GR_SUCCESS && gr_mpoly_equal(q, f, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: quotient != f\n");
                        flint_printf("alg = %d\n", alg);
                        gr_ctx_println(ctx);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(h, ctx);
        gr_mpoly_clear(q, ctx);
        gr_mpoly_clear(r, ctx);
        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
        gr_ctx_clear(ctx1);
    }

    TEST_FUNCTION_END(state);
}
