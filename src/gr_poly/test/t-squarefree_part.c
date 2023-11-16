/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_squarefree_part, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx, poly_ctx, fmpz_ctx;
        gr_poly_t A, B, C;
        gr_ptr c;
        gr_vec_t F, exp;
        slong i;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        while (gr_ctx_is_field(ctx) != T_TRUE)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_ctx_init_gr_poly(poly_ctx, ctx);
        gr_ctx_init_fmpz(fmpz_ctx);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);

        gr_vec_init(F, 0, poly_ctx);
        gr_vec_init(exp, 0, fmpz_ctx);

        c = gr_heap_init(ctx);

        if (ctx->methods == _ca_methods)
            status |= gr_poly_randtest(A, state, 3, ctx);
        else
            status |= gr_poly_randtest(A, state, 10, ctx);

        status |= gr_poly_factor_squarefree(c, F, exp, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_squarefree_part(B, A, ctx);

            if (status == GR_SUCCESS)
            {
                status |= gr_poly_one(C, ctx);
                for (i = 0; i < F->length; i++)
                    status |= gr_poly_mul(C, C, gr_vec_entry_ptr(F, i, poly_ctx), ctx);

                if (status == GR_SUCCESS && gr_poly_equal(B, C, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (product)\n\n");
                    flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                    flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
                    flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");

                    flint_abort();
                }
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(C, ctx);

        gr_vec_clear(F, poly_ctx);
        gr_vec_clear(exp, fmpz_ctx);

        gr_heap_clear(c, ctx);

        gr_ctx_clear(poly_ctx);
        gr_ctx_clear(fmpz_ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
