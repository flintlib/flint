/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

FLINT_DLL extern gr_static_method_table _ca_methods;

void _gr_poly_test_gcd_field(gr_method_poly_gcd_op gcd_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        if (gr_ctx_is_field(ctx) == T_TRUE)
        {
            gr_poly_t A, B, C, AC, BC, MC, G;
            int status = GR_SUCCESS;
            slong n;

            gr_poly_init(A, ctx);
            gr_poly_init(B, ctx);
            gr_poly_init(C, ctx);
            gr_poly_init(AC, ctx);
            gr_poly_init(BC, ctx);
            gr_poly_init(MC, ctx);
            gr_poly_init(G, ctx);

            n = 1 + n_randint(state, maxn);

            /* Hack: slow */
            if (ctx->methods == _ca_methods)
                n = 3;

            status = gr_poly_randtest(A, state, n, ctx);
            status |= gr_poly_randtest(B, state, n, ctx);

            status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, A, B, ctx);

            if (status == GR_SUCCESS)
            {
                if (gr_poly_is_zero(G, ctx) == T_FALSE)
                {
                    status |= gr_poly_divrem(AC, BC, A, G, ctx);

                    if (status == GR_SUCCESS && gr_poly_is_zero(BC, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: gcd does not divide A\n\n");
                        gr_ctx_println(ctx);
                        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_abort();
                    }

                    status |= gr_poly_divrem(AC, BC, B, G, ctx);

                    if (status == GR_SUCCESS && gr_poly_is_zero(BC, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: gcd does not divide B\n\n");
                        gr_ctx_println(ctx);
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
                        status |= gr_poly_gcd(G, G, BC, ctx);
                        break;
                    case 1:
                        status |= gr_poly_set(G, BC, ctx);
                        status |= gr_poly_gcd(G, AC, G, ctx);
                        break;
                    default:
                        status |= gr_poly_gcd(G, AC, BC, ctx);
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
                            gr_ctx_println(ctx);
                            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                            flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                            flint_abort();
                        }
                    }
                }
            }

            status = gr_poly_randtest(A, state, n, ctx);
            status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, A, A, ctx);

            if (status == GR_SUCCESS)
            {
                status |= gr_poly_make_monic(B, A, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(G, B, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (self)\n\n");
                    gr_ctx_println(ctx);
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
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}

void _gr_poly_test_gcd_ufd(gr_method_poly_gcd_op gcd_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        if (gr_ctx_is_unique_factorization_domain(ctx) == T_TRUE)
        {
            gr_poly_t A, B, C, AC, BC, MC, G, R;
            int status = GR_SUCCESS;
            slong n;

            gr_poly_init(A, ctx);
            gr_poly_init(B, ctx);
            gr_poly_init(C, ctx);
            gr_poly_init(AC, ctx);
            gr_poly_init(BC, ctx);
            gr_poly_init(MC, ctx);
            gr_poly_init(G, ctx);
            gr_poly_init(R, ctx);

            n = 1 + n_randint(state, maxn);
            if (ctx->methods == _ca_methods)
                n = 3;

            status = gr_poly_randtest(A, state, n, ctx);
            status |= gr_poly_randtest(B, state, n, ctx);
            status |= gr_poly_randtest(C, state, n, ctx);

            status |= gr_poly_mul(AC, A, C, ctx);
            status |= gr_poly_mul(BC, B, C, ctx);

            switch (n_randint(state, 3))
            {
                case 0:
                    status |= gr_poly_set(G, AC, ctx);
                    status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, G, BC, ctx);
                    break;
                case 1:
                    status |= gr_poly_set(G, BC, ctx);
                    status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, AC, G, ctx);
                    break;
                default:
                    status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, AC, BC, ctx);
                    break;
            }

            /* Check that C divides G = gcd(AC, BC) */

            if (status == GR_SUCCESS)
            {
                if (gr_poly_is_zero(C, ctx) == T_FALSE)
                {
                    status |= gr_poly_rem(R, G, C, ctx);

                    if (status == GR_SUCCESS && gr_poly_is_zero(R, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: C does not divide gcd(AC, BC)\n\n");
                        gr_ctx_println(ctx);
                        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                        flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                        flint_printf("AC = "); gr_poly_print(AC, ctx); flint_printf("\n");
                        flint_printf("BC = "); gr_poly_print(BC, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_printf("R = "); gr_poly_print(R, ctx); flint_printf("\n");
                        flint_abort();
                    }
                }

                if (gr_poly_is_zero(G, ctx) == T_FALSE)
                {
                    status |= gr_poly_rem(R, AC, G, ctx);

                    if (status == GR_SUCCESS && gr_poly_is_zero(R, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: gcd does not divide A\n\n");
                        gr_ctx_println(ctx);
                        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                        flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                        flint_printf("AC = "); gr_poly_print(AC, ctx); flint_printf("\n");
                        flint_printf("BC = "); gr_poly_print(BC, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_printf("R = "); gr_poly_print(R, ctx); flint_printf("\n");
                        flint_abort();
                    }

                    status |= gr_poly_rem(R, BC, G, ctx);

                    if (status == GR_SUCCESS && gr_poly_is_zero(R, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL: gcd does not divide BC\n\n");
                        gr_ctx_println(ctx);
                        flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                        flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                        flint_printf("AC = "); gr_poly_print(AC, ctx); flint_printf("\n");
                        flint_printf("BC = "); gr_poly_print(BC, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_printf("R = "); gr_poly_print(R, ctx); flint_printf("\n");
                        flint_abort();
                    }
                }
            }

            /* gcd(A,B) = 1  ==>  gcd(AC, BC) = canonical_associate(C) */
            status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, A, B, ctx);

            if (status == GR_SUCCESS && gr_poly_is_one(G, ctx) == T_TRUE)
            {
                status |= gr_poly_randtest(C, state, n, ctx);

                status |= gr_poly_mul(AC, A, C, ctx);
                status |= gr_poly_mul(BC, B, C, ctx);

                switch (n_randint(state, 3))
                {
                    case 0:
                        status |= gr_poly_set(G, AC, ctx);
                        status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, G, BC, ctx);
                        break;
                    case 1:
                        status |= gr_poly_set(G, BC, ctx);
                        status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, AC, G, ctx);
                        break;
                    default:
                        status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, AC, BC, ctx);
                        break;
                }

                if (status == GR_SUCCESS)
                {
                    if (gr_poly_is_zero(C, ctx) == T_FALSE)
                    {
                        status |= gr_poly_canonical_associate(MC, NULL, C, ctx);

                        if (status == GR_SUCCESS && gr_poly_equal(MC, G, ctx) == T_FALSE)
                        {
                            flint_printf("FAIL\n\n");
                            gr_ctx_println(ctx);
                            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                            flint_printf("MC = "); gr_poly_print(MC, ctx); flint_printf("\n");
                            flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                            flint_abort();
                        }
                    }
                }
            }

            if ((ctx->which_ring == GR_CTX_FMPQ ||
                 ctx->which_ring == GR_CTX_FMPZ ||
                 ctx->which_ring == GR_CTX_FMPZI ||
                 ctx->which_ring == GR_CTX_FMPZ_POLY ||
                 ctx->which_ring == GR_CTX_FMPQ_POLY ||
                 ctx->which_ring == GR_CTX_FMPZ_MPOLY ||
                (gr_ctx_is_field(ctx) == T_TRUE && gr_ctx_is_finite(ctx) == T_TRUE)) && status != GR_SUCCESS)
            {
                flint_printf("FAIL: unexpectedly did not succeed\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
                flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                flint_abort();
            }

            status = gr_poly_randtest(A, state, n, ctx);
            status |= gr_poly_gcd_wrapper(gcd_impl, 1, G, A, A, ctx);

            if (status == GR_SUCCESS)
            {
                status |= gr_poly_canonical_associate(B, NULL, A, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(G, B, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (self gcd)\n\n");
                    gr_ctx_println(ctx);
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
            gr_poly_clear(R, ctx);
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}

POP_OPTIONS
