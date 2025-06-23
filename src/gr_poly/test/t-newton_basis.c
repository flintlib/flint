/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_newton_basis, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_vec_t R, Y;
        gr_poly_t F, G, H, A, B;
        gr_ptr x, y, z;
        slong i;

        gr_ctx_init_random(ctx, state);
        /* Hack: avoid because slow */
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_vec_init(R, n_randint(state, 6), ctx);
        gr_poly_init(F, ctx);
        gr_poly_init(G, ctx);
        gr_poly_init(H, ctx);
        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);

        status |= _gr_vec_randtest(R->entries, state, R->length, ctx);
        status |= gr_poly_randtest(F, state, n_randint(state, 6), ctx);
        status |= gr_poly_randtest(G, state, n_randint(state, 6), ctx);
        status |= gr_poly_randtest(H, state, n_randint(state, 6), ctx);
        status |= gr_randtest(x, state, ctx);
        status |= gr_randtest(y, state, ctx);
        status |= gr_randtest(z, state, ctx);

        status |= gr_poly_newton_basis_from_monomial(G, R, F, ctx);

        if (status == GR_SUCCESS)
        {
            status = gr_poly_evaluate(y, F, x, ctx);
            status |= gr_poly_newton_basis_evaluate(z, R, G, x, ctx);

            if (status == GR_SUCCESS && gr_equal(y, z, ctx) == T_FALSE)
            {
                flint_printf("FAIL (evaluation)\n\n");
                gr_ctx_println(ctx);
                flint_printf("R = "); gr_vec_print(R, ctx); flint_printf("\n");
                flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                flint_printf("x = %{gr}\n\n", x);
                flint_printf("F(x) = %{gr}\n\n", y);
                flint_printf("G(x) = %{gr}\n\n", z);
                flint_abort();
            }

            status = gr_poly_newton_basis_to_monomial(H, R, G, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(F, H, ctx) == T_FALSE)
            {
                flint_printf("FAIL (roundtrip)\n\n");
                gr_ctx_println(ctx);
                flint_printf("R = "); gr_vec_print(R, ctx); flint_printf("\n");
                flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                flint_printf("H = "); gr_poly_print(G, ctx); flint_printf("\n");
                flint_abort();
            }

            if (R->length >= F->length && gr_ctx_is_integral_domain(ctx) == T_TRUE)
            {
                gr_vec_init(Y, F->length, ctx);
                for (i = 0; i < F->length; i++)
                    status |= gr_poly_newton_basis_evaluate(gr_vec_entry_ptr(Y, i, ctx),
                                        R, G, gr_vec_entry_srcptr(R, i, ctx), ctx);

                status |= gr_poly_randtest(A, state, n_randint(state, 5), ctx);
                status |= gr_poly_randtest(B, state, n_randint(state, 5), ctx);

                status |= gr_poly_newton_basis_interpolate(A, R, Y, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(A, G, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (interpolation)\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("R = "); gr_vec_print(R, ctx); flint_printf("\n");
                    flint_printf("Y = "); gr_vec_print(Y, ctx); flint_printf("\n");
                    flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                    flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                    flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                    flint_abort();
                }

                if (status == GR_SUCCESS)
                {
                    status |= gr_poly_newton_basis_interpolate_exact(B, R, Y, ctx);

                    if (status != GR_SUCCESS || gr_poly_equal(B, G, ctx) == T_FALSE)
                    {
                        flint_printf("FAIL (interpolation 2)\n\n");
                        gr_ctx_println(ctx);
                        flint_printf("R = "); gr_vec_print(R, ctx); flint_printf("\n");
                        flint_printf("Y = "); gr_vec_print(Y, ctx); flint_printf("\n");
                        flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                        flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                        flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                        flint_abort();
                    }
                }

                gr_vec_clear(Y, ctx);
            }
        }
        else if (R->length >= F->length - 1 && ctx->which_ring == GR_CTX_FMPZ)
        {
            flint_printf("FAIL (expected success)\n\n");
            gr_ctx_println(ctx);
            flint_printf("R = "); gr_vec_print(R, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_vec_clear(R, ctx);
        gr_poly_clear(F, ctx);
        gr_poly_clear(G, ctx);
        gr_poly_clear(H, ctx);
        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
