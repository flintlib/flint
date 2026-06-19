/*
    Copyright (C) 2023, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_squarefree_part, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t A, B, C;
        gr_ptr c;
        gr_poly_vec_t F;
        fmpz_vec_t exp;
        slong i;
        int status = GR_SUCCESS;

        gr_ctx_init_random_commutative_ring(ctx, state);

        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random_commutative_ring(ctx, state);
        }

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(C, ctx);

        gr_poly_vec_init(F, 0, ctx);
        fmpz_vec_init(exp, 0);

        c = gr_heap_init(ctx);

        status |= gr_poly_randtest(A, state, 8, ctx);

        status |= gr_poly_factor_squarefree(c, F, exp, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_squarefree_part(B, A, ctx);

            if (status == GR_SUCCESS)
            {
                status |= gr_poly_one(C, ctx);
                for (i = 0; i < F->length; i++)
                    status |= gr_poly_mul(C, C, F->entries + i, ctx);
                status |= gr_poly_canonical_associate(C, NULL, C, ctx);

                if (status == GR_SUCCESS && gr_poly_equal(B, C, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (product)\n\n");
                    gr_ctx_println(ctx);
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

        gr_poly_vec_clear(F, ctx);
        fmpz_vec_clear(exp);
        gr_heap_clear(c, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
