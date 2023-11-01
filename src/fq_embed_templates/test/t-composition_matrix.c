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

TEST_TEMPLATE_FUNCTION_START(T, embed_composition_matrix, state)
{
    int i;

    /* Check that Mat(a^p) = Mat(x^p) * Mat(a) for random a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) frob, a;
        const TEMPLATE(B, poly_struct) *modulus;
        TEMPLATE(B, mat_t) mat_frob, mat_a, mat_aq, res;
        slong d;

        while (TEMPLATE(T, ctx_randtest)(ctx, state),
               d = TEMPLATE(T, ctx_degree)(ctx),
               d == 1)
        {
            TEMPLATE(T, ctx_clear)(ctx);
        }

        modulus = TEMPLATE(T, ctx_modulus)(ctx);

        TEMPLATE(T, init)(frob, ctx);
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(B, mat_init)(mat_frob, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(mat_a, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(mat_aq, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(res, d, d, TEMPLATE(B, poly_modulus)(modulus));

        TEMPLATE(T, gen)(frob, ctx);
        TEMPLATE(T, pow)(frob, frob, TEMPLATE(T, ctx_prime)(ctx), ctx);
        TEMPLATE(T, embed_composition_matrix)(mat_frob, frob, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, embed_composition_matrix)(mat_a, a, ctx);

        TEMPLATE(B, mat_mul)(res, mat_frob, mat_a);

        TEMPLATE(T, pow)(a, a, TEMPLATE(T, ctx_prime)(ctx), ctx);
        TEMPLATE(T, embed_composition_matrix)(mat_aq, a, ctx);

        if (!TEMPLATE(B, mat_equal)(res, mat_aq))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("x^q: "), TEMPLATE(T, print_pretty)(frob, ctx), flint_printf("\n");
            flint_printf("M(x^q)*M(a) = M(a^q)\n"),
                TEMPLATE(B, mat_print_pretty)(mat_frob), flint_printf("\n"),
                TEMPLATE(B, mat_print_pretty)(mat_a), flint_printf("\n"),
                TEMPLATE(B, mat_print_pretty)(mat_aq), flint_printf("\n"),
                TEMPLATE(B, mat_print_pretty)(res), flint_printf("\n");

            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(B, mat_clear)(mat_frob);
        TEMPLATE(B, mat_clear)(mat_a);
        TEMPLATE(B, mat_clear)(mat_aq);
        TEMPLATE(B, mat_clear)(res);
        TEMPLATE(T, clear)(frob, ctx);
        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
#endif
