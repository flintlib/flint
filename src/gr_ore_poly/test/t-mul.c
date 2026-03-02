/*
    Copyright (C) 2026 Maria Neagoie, supervised by Marc Mezzarobba and Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_ore_poly.h"
#include "gr_vec.h"

// same as GR_MUST_SUCCEED but details of the failing printed
#define MUST_OK(expr) do { \
    int st = (expr); \
    if (st != GR_SUCCESS) \
    { \
        flint_printf("FAIL at %s:%d: %s returned %d\n", __FILE__, __LINE__, #expr, st); \
        flint_abort(); \
    } \
} while(0)

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_ore_poly_mul, state)
{
    for (slong iter = 0; iter < 1500 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_ore_poly_ctx_t ore_ctx;
        slong reps;
        ore_algebra_t algebra;
        int found_algebra = 0;

        while(!found_algebra)
        {
            // Random context
            gr_ore_poly_ctx_init_randtest2(ctx, ore_ctx, state);

            // Random Algebra
            algebra = GR_ORE_POLY_CTX(ore_ctx)->which_algebra;

            if (algebra == ORE_ALGEBRA_COMMUTATIVE) // Standard polynomials
                found_algebra = 1;
            else
            {
                if (algebra != ORE_ALGEBRA_FROBENIUS && ctx->which_ring == GR_CTX_GR_POLY)
                {
                    gr_ctx_struct* base_ring = POLYNOMIAL_CTX(ctx)->base_ring;
                    gr_ctx_clear(ctx);
                    gr_ctx_init_gr_poly(ctx, base_ring);
                    found_algebra = 1;
                }
            }
        }

        // gr_ctx_print(ctx);
        // flint_printf("\n");

        // flint_printf("sigma_delta = %p\n", (void*) GR_ORE_POLY_CTX(ore_ctx)->sigma_delta);
        /*if (GR_ORE_POLY_CTX(ore_ctx)->sigma_delta == NULL)
        {
            gr_ore_poly_ctx_clear(ore_ctx);
            gr_ctx_clear(ctx);
            continue;
        }*/

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

        /* Hack: for string conversion tests, make sure we don't have
           overlapping generator names. */
        gr_vec_t vec;
        gr_vec_init(vec, 0, ctx);
        if (gr_gens_recursive(vec, ctx) == GR_SUCCESS)
        {
            const char * vars[] = { "DD" };

            MUST_OK(gr_ctx_set_gen_names(ore_ctx, vars));

        }
        gr_vec_clear(vec, ctx);

        for (slong j = 0; j < reps; j++)
        {
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
            MUST_OK(gr_ore_poly_randtest(A, state, 4, ore_ctx));
            MUST_OK(gr_ore_poly_randtest(B, state, 4, ore_ctx));
            MUST_OK(gr_ore_poly_randtest(C, state, 4, ore_ctx));

            // flint_printf("A = "); gr_ore_poly_print(A, ore_ctx); flint_printf("\n");
            // flint_printf("B = "); gr_ore_poly_print(B, ore_ctx); flint_printf("\n");
            // flint_printf("C = "); gr_ore_poly_print(C, ore_ctx); flint_printf("\n");
            // flint_printf("\n");

            MUST_OK(gr_ore_poly_mul(AB, A, B, ore_ctx));
            // Aliasing test
            if (iter & 1) // odd => test out = A
            {
                MUST_OK(gr_ore_poly_set(TEMP, A, ore_ctx));
                MUST_OK(gr_ore_poly_mul(TEMP, TEMP, B, ore_ctx));
            }
            else // even => test out = B
            {
                MUST_OK(gr_ore_poly_set(TEMP, B, ore_ctx));
                MUST_OK(gr_ore_poly_mul(TEMP, A, TEMP, ore_ctx));
            }

            if (gr_ore_poly_equal(TEMP, AB, ore_ctx) == T_FALSE)
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
            MUST_OK(gr_ore_poly_mul(LHS, AB, C, ore_ctx));

            MUST_OK(gr_ore_poly_mul(BC, B, C, ore_ctx));
            MUST_OK(gr_ore_poly_mul(RHS, A, BC, ore_ctx));

            if (gr_ore_poly_equal(LHS, RHS, ore_ctx) == T_FALSE)
            {
                flint_printf("FAIL: (A·B)·C = A·(B·C)\n");

                flint_printf("A = "); gr_ore_poly_print(A, ore_ctx); flint_printf("\n");
                flint_printf("B = "); gr_ore_poly_print(B, ore_ctx); flint_printf("\n");
                flint_printf("C = "); gr_ore_poly_print(C, ore_ctx); flint_printf("\n");

                flint_printf("(A·B)·C = "); gr_ore_poly_print(LHS, ore_ctx); flint_printf("\n");
                flint_printf("A·(B·C) = "); gr_ore_poly_print(RHS, ore_ctx); flint_printf("\n");

                flint_abort();
            }

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

    TEST_FUNCTION_END(state);
}
