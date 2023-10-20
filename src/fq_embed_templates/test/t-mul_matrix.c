/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, embed_mul_matrix, state)
{
    int i;

    /* Check that Mat(a^2) = Mat(a)^2 for random a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a;
        TEMPLATE(B, mat_t) mat_a, mat_a_sq, mat_a_a;
        const TEMPLATE(B, poly_struct) *modulus;
        slong d;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        d = TEMPLATE(T, ctx_degree)(ctx);
        modulus = TEMPLATE(T, ctx_modulus)(ctx);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(B, mat_init)(mat_a, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(mat_a_sq, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(mat_a_a, d, d, TEMPLATE(B, poly_modulus)(modulus));

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, embed_mul_matrix)(mat_a, a, ctx);

        TEMPLATE(T, mul)(a, a, a, ctx);
        TEMPLATE(T, embed_mul_matrix)(mat_a_sq, a, ctx);

        TEMPLATE(B, mat_mul)(mat_a_a, mat_a, mat_a);

        if (!TEMPLATE(B, mat_equal)(mat_a_a, mat_a_sq)) {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("a^2: "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("M(a)^2 = M(a^2)\n"),
                TEMPLATE(B, mat_print_pretty)(mat_a), flint_printf("^2\n=\n"),
                TEMPLATE(B, mat_print_pretty)(mat_a_a), flint_printf("\n=\n"),
                TEMPLATE(B, mat_print_pretty)(mat_a_sq), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(B, mat_clear)(mat_a);
        TEMPLATE(B, mat_clear)(mat_a_sq);
        TEMPLATE(B, mat_clear)(mat_a_a);
        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
#endif
