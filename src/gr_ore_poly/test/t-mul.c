/*
    Copyright (C) 2026 Maria Neagoie

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_ore_poly.h"
#include "gr_vec.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_ore_poly_mul, state)
{
    int status;
    // slong success = 0, unable = 0, unable_frobenius = 0, domain = 0, success_frobenius = 0;
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_ore_poly_ctx_t ore_ctx;
        slong reps;
        ore_algebra_t algebra;

        while(1)
        {
            // Random context
            gr_ore_poly_ctx_init_randtest2(ctx, ore_ctx, state);

            // Random Algebra
            algebra = GR_ORE_POLY_CTX(ore_ctx)->which_algebra;

            if (algebra == ORE_ALGEBRA_COMMUTATIVE) // Standard polynomials
                break;

            if (algebra == ORE_ALGEBRA_FROBENIUS && gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
                break;

            if (algebra != ORE_ALGEBRA_FROBENIUS && ctx->which_ring == GR_CTX_GR_POLY)
                break;

            gr_ore_poly_ctx_clear(ore_ctx);
            gr_ctx_clear(ctx);
        }

        if (algebra == ORE_ALGEBRA_MAHLER || algebra == ORE_ALGEBRA_Q_SHIFT)
            ore_ctx->size_limit = 4;

        if (gr_ctx_is_finite(ctx) == T_TRUE ||
            gr_ctx_has_real_prec(ctx) == T_TRUE)
        {
            reps = 10;
        }
        else if (ctx->methods == _ca_methods) /* hack: slow */
        {
            reps = 1;
        }
        else
        {
            reps = 3;
        }

        for (slong j = 0; j < reps; j++)
        {
            status = GR_SUCCESS;

            gr_ore_poly_t A, B, C, AB, BC, LHS, RHS, TEMP;

            gr_ore_poly_init(A, ore_ctx);
            gr_ore_poly_init(B, ore_ctx);
            gr_ore_poly_init(C, ore_ctx);
            gr_ore_poly_init(AB, ore_ctx);
            gr_ore_poly_init(BC, ore_ctx);
            gr_ore_poly_init(LHS, ore_ctx);
            gr_ore_poly_init(RHS, ore_ctx);
            gr_ore_poly_init(TEMP, ore_ctx);

            // Build objects with coefficients in random context
            status |= gr_ore_poly_randtest(A, state, 1 + n_randint(state, 6), ore_ctx);
            status |= gr_ore_poly_randtest(B, state, 1 + n_randint(state, 6), ore_ctx);
            status |= gr_ore_poly_randtest(C, state, 1 + n_randint(state, 6), ore_ctx);

            status |= gr_ore_poly_mul(AB, A, B, ore_ctx);
            // Aliasing test
            if (iter & 1) // odd => test out = A
            {
                status |= gr_ore_poly_set(TEMP, A, ore_ctx);
                status |= gr_ore_poly_mul(TEMP, TEMP, B, ore_ctx);
            }
            else // even => test out = B
            {
                status |= gr_ore_poly_set(TEMP, B, ore_ctx);
                status |= gr_ore_poly_mul(TEMP, A, TEMP, ore_ctx);
            }

            if (status == GR_SUCCESS && gr_ore_poly_equal(TEMP, AB, ore_ctx) == T_FALSE)
            {
                if(iter & 1)
                    flint_printf("FAIL: alias (out = A)\n");
                else
                    flint_printf("FAIL: alias (out = B)\n");

                flint_printf("A = "); gr_ore_poly_print(A, ore_ctx); flint_printf("\n");
                flint_printf("B = "); gr_ore_poly_print(B, ore_ctx); flint_printf("\n");

                flint_printf("A·B = "); gr_ore_poly_print(AB, ore_ctx); flint_printf("\n");
                flint_printf("TEMP = "); gr_ore_poly_print(TEMP, ore_ctx); flint_printf("\n");

                flint_abort();
            }

            // Associativity test: (A*B)*C = A*(B*C)
            status |= gr_ore_poly_mul(LHS, AB, C, ore_ctx);

            status |= gr_ore_poly_mul(BC, B, C, ore_ctx);
            status |= gr_ore_poly_mul(RHS, A, BC, ore_ctx);

            if (status == GR_SUCCESS && gr_ore_poly_equal(LHS, RHS, ore_ctx) == T_FALSE)
            {
                flint_printf("FAIL: (A·B)·C = A·(B·C)\n");

                flint_printf("A = "); gr_ore_poly_print(A, ore_ctx); flint_printf("\n");
                flint_printf("B = "); gr_ore_poly_print(B, ore_ctx); flint_printf("\n");
                flint_printf("C = "); gr_ore_poly_print(C, ore_ctx); flint_printf("\n");

                flint_printf("(A·B)·C = "); gr_ore_poly_print(LHS, ore_ctx); flint_printf("\n");
                flint_printf("A·(B·C) = "); gr_ore_poly_print(RHS, ore_ctx); flint_printf("\n");

                flint_abort();
            }

            /*if(status == GR_SUCCESS)
            {
                success++;
                if (algebra == ORE_ALGEBRA_FROBENIUS)
                    success_frobenius++;
            }

            if (status == GR_DOMAIN)
                domain++;

            if (status == GR_UNABLE)
            {
                unable++;
                if (algebra == ORE_ALGEBRA_FROBENIUS)
                {
                    unable_frobenius++;
                    gr_ctx_println(ctx);
                }
            }*/

            gr_ore_poly_clear(A, ore_ctx);
            gr_ore_poly_clear(B, ore_ctx);
            gr_ore_poly_clear(C, ore_ctx);
            gr_ore_poly_clear(AB, ore_ctx);
            gr_ore_poly_clear(BC, ore_ctx);
            gr_ore_poly_clear(LHS, ore_ctx);
            gr_ore_poly_clear(RHS, ore_ctx);
            gr_ore_poly_clear(TEMP, ore_ctx);
        }

        gr_ore_poly_ctx_clear(ore_ctx);
        gr_ctx_clear(ctx);
    }
    /*flint_printf("GR_SUCCESS = %d\n", success);
    flint_printf("GR_DOMAIN = %d\n", domain);
    flint_printf("GR_UNABLE = %d\n", unable);
    flint_printf("GR_UNABLE Frobenius = %d\n", unable_frobenius);
    flint_printf("GR_SUCCESS Frobenius = %d\n", success_frobenius);*/
    TEST_FUNCTION_END(state);
}
