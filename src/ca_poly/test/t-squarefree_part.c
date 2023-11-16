/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca_poly.h"

TEST_FUNCTION_START(ca_poly_squarefree_part, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_poly_t A, B, C;
        ca_t c;
        ca_poly_vec_t F;
        ulong * exp;
        slong i;

        ca_ctx_init(ctx);

        ca_poly_init(A, ctx);
        ca_poly_init(B, ctx);
        ca_poly_init(C, ctx);
        ca_poly_vec_init(F, 0, ctx);
        ca_init(c, ctx);

        ca_poly_randtest(A, state, 5, 1, 5, ctx);
        exp = flint_malloc(sizeof(ulong) * A->length);

        if (ca_poly_factor_squarefree(c, F, exp, A, ctx))
        {
            if (ca_poly_squarefree_part(B, A, ctx))
            {
                ca_poly_one(C, ctx);
                for (i = 0; i < F->length; i++)
                        ca_poly_mul(C, C, F->entries + i, ctx);

                if (ca_poly_check_equal(B, C, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (product)\n\n");
                    flint_printf("A = "); ca_poly_print(A, ctx); flint_printf("\n");
                    flint_printf("B = "); ca_poly_print(B, ctx); flint_printf("\n");
                    flint_printf("C = "); ca_poly_print(C, ctx); flint_printf("\n");

                    flint_abort();
                }
            }
        }

        flint_free(exp);

        ca_poly_clear(A, ctx);
        ca_poly_clear(B, ctx);
        ca_poly_clear(C, ctx);
        ca_poly_vec_clear(F, ctx);
        ca_clear(c, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
