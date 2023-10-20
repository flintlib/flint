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

TEST_TEMPLATE_FUNCTION_START(T, embed_mono_dual_matrix, state)
{
    int i;

    /* Check that the two functions are inverse of one another */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        const TEMPLATE(B, poly_struct) *modulus;
        TEMPLATE(B, mat_t) m2d, d2m, one, two;
        slong d;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        d = TEMPLATE(T, ctx_degree)(ctx);
        modulus = TEMPLATE(T, ctx_modulus)(ctx);

        TEMPLATE(B, mat_init)(m2d, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(d2m, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(one, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(two, d, d, TEMPLATE(B, poly_modulus)(modulus));

        TEMPLATE(T, embed_mono_to_dual_matrix)(m2d, ctx);
        TEMPLATE(T, embed_dual_to_mono_matrix)(d2m, ctx);
        TEMPLATE(B, mat_mul)(one, m2d, d2m);

        TEMPLATE(B, mat_one)(two);

        if (!TEMPLATE(B, mat_equal)(one, two)) {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("Mono -> Dual\n"),
                TEMPLATE(B, mat_print_pretty)(m2d), flint_printf("\nDual -> Mono\n"),
                TEMPLATE(B, mat_print_pretty)(d2m), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(B, mat_clear)(m2d);
        TEMPLATE(B, mat_clear)(d2m);
        TEMPLATE(B, mat_clear)(one);
        TEMPLATE(B, mat_clear)(two);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
#endif
