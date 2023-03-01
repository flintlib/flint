/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("gcd_euclidean....");
    fflush(stdout);

    flint_randinit(state);


    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, C, AC, BC, MC, G;
        int status = GR_SUCCESS;
        slong n;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);
        gr_poly_init(AC, ctx);
        gr_poly_init(BC, ctx);
        gr_poly_init(MC, ctx);
        gr_poly_init(G, ctx);

        n = 6;
        if (ctx->which_ring == GR_CTX_CC_CA || ctx->which_ring == GR_CTX_RR_CA)
            n = 3;

        status = gr_poly_randtest(A, state, n, ctx);
        status |= gr_poly_randtest(B, state, n, ctx);
        status |= gr_poly_gcd_euclidean(G, A, B, ctx);

        if (status == GR_SUCCESS)
        {
            if (gr_poly_is_zero(G, ctx) == T_FALSE)
            {
                status |= gr_poly_divrem(AC, BC, A, G, ctx);

                if (status == GR_SUCCESS && gr_poly_is_zero(BC, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: gcd does not divide A\n\n");
                    flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                    flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                    flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                    flint_abort();
                }

                status |= gr_poly_divrem(AC, BC, B, G, ctx);

                if (status == GR_SUCCESS && gr_poly_is_zero(BC, ctx) == T_FALSE)
                {
                    flint_printf("FAIL: gcd does not divide B\n\n");
                    flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                    flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                    flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        if (status == GR_SUCCESS && gr_poly_is_one(G, ctx) == T_TRUE)
        {
            status |= gr_poly_randtest(C, state, n, ctx);

            status |= gr_poly_mul(AC, A, C, ctx);
            status |= gr_poly_mul(BC, B, C, ctx);

            switch (n_randint(state, 3))
            {
                case 0:
                    status |= gr_poly_set(G, AC, ctx);
                    status |= gr_poly_gcd_euclidean(G, G, BC, ctx);
                    break;
                case 1:
                    status |= gr_poly_set(G, BC, ctx);
                    status |= gr_poly_gcd_euclidean(G, AC, G, ctx);
                    break;
                default:
                    status |= gr_poly_gcd_euclidean(G, AC, BC, ctx);
                    break;
            }

            if (status == GR_SUCCESS)
            {
                if (gr_poly_is_zero(C, ctx) == T_FALSE)
                {
                    status |= gr_poly_make_monic(MC, C, ctx);

                    if (status == GR_SUCCESS && gr_poly_equal(MC, G, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL\n\n");
                        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                        flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_abort();
                    }
                }
            }
        }

        if ((ctx->which_ring == GR_CTX_FMPQ || (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE)) && status != GR_SUCCESS)
        {
            flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
            flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
            flint_abort();
        }

        status = gr_poly_randtest(A, state, n, ctx);
        status |= gr_poly_gcd_euclidean(G, A, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_make_monic(B, A, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(G, B, ctx) == T_FALSE)
            {
                flint_printf("FAIL (self)\n\n");
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(C, ctx);
        gr_poly_clear(AC, ctx);
        gr_poly_clear(BC, ctx);
        gr_poly_clear(MC, ctx);
        gr_poly_clear(G, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
