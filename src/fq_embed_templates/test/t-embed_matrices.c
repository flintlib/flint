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

TEST_TEMPLATE_FUNCTION_START(T, embed_matrices, state)
{
    int i, j;

    /* Check that isomorphism to self gives identity matrices */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) gen;
        const TEMPLATE(B, poly_struct) *modulus;
        TEMPLATE(B, mat_t) embed, project, one;
        slong d;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        d = TEMPLATE(T, ctx_degree)(ctx);
        modulus = TEMPLATE(T, ctx_modulus)(ctx);

        TEMPLATE(T, init)(gen, ctx);
        TEMPLATE(T, gen)(gen, ctx);
        TEMPLATE(T, pow)(gen, gen, TEMPLATE(T, ctx_prime)(ctx), ctx);

        TEMPLATE(B, mat_init)(embed, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(project, d, d, TEMPLATE(B, poly_modulus)(modulus));
        TEMPLATE(B, mat_init)(one, d, d, TEMPLATE(B, poly_modulus)(modulus));

        TEMPLATE(T, embed_matrices)(embed, project, gen, ctx, gen, ctx, modulus);
        TEMPLATE(B, mat_one)(one);

        if (!TEMPLATE(B, mat_equal)(embed, one) || !TEMPLATE(B, mat_equal)(project, one)) {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("x^p: "), TEMPLATE(T, print_pretty)(gen, ctx), flint_printf("\n");
            flint_printf("Embed\n"),
                TEMPLATE(B, mat_print_pretty)(embed), flint_printf("\nProject\n"),
                TEMPLATE(B, mat_print_pretty)(project), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(B, mat_clear)(embed);
        TEMPLATE(B, mat_clear)(project);
        TEMPLATE(B, mat_clear)(one);
        TEMPLATE(T, clear)(gen, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check random emebedding (degrees 1..5) */
    for (j = 1; j < 6; j++) {
        for (i = 0; i < (6 - j) * flint_test_multiplier(); i++)
        {
            TEMPLATE(T, ctx_t) ctx1, ctx2;
            TEMPLATE(T, t) gen1, gen2;
            TEMPLATE(B, poly_t) minpoly;
            const TEMPLATE(B, poly_struct) *modulus;
            TEMPLATE(B, poly_t) modulus2;
            TEMPLATE(B, mat_t) embed, project, comp, one;
            slong m, n;

            while (TEMPLATE(T, ctx_randtest)(ctx1, state),
                    m = TEMPLATE(T, ctx_degree)(ctx1),
                    m == 1)
            {
                TEMPLATE(T, ctx_clear)(ctx1);
            }

            n = m*j;

            modulus = TEMPLATE(T, ctx_modulus)(ctx1);

            TEMPLATE(B, poly_init)(modulus2, TEMPLATE(B, poly_modulus)(modulus));
            TEMPLATE(B, poly_randtest_monic_irreducible)(modulus2, state, n+1);
            TEMPLATE(T, ctx_init_modulus)(ctx2, modulus2, "X");

            TEMPLATE(T, init)(gen1, ctx1);
            TEMPLATE(T, init)(gen2, ctx2);
            TEMPLATE(B, poly_init)(minpoly, TEMPLATE(B, poly_modulus)(modulus));
            TEMPLATE(T, embed_gens)(gen1, gen2, minpoly, ctx1, ctx2);

            TEMPLATE(B, mat_init)(embed, n, m, TEMPLATE(B, poly_modulus)(modulus));
            TEMPLATE(B, mat_init)(project, m, n, TEMPLATE(B, poly_modulus)(modulus));
            TEMPLATE(B, mat_init)(comp, m, m, TEMPLATE(B, poly_modulus)(modulus));
            TEMPLATE(B, mat_init)(one, m, m, TEMPLATE(B, poly_modulus)(modulus));

            TEMPLATE(T, embed_matrices)(embed, project, gen1, ctx1, gen2, ctx2, minpoly);

            TEMPLATE(B, mat_mul)(comp, project, embed);
            TEMPLATE(B, mat_one)(one);
            if (!TEMPLATE(B, mat_equal)(comp, one)) {
                flint_printf("FAIL:\n\n");
                flint_printf("CTX 1\n"), TEMPLATE(T, ctx_print)(ctx1), flint_printf("\n");
                flint_printf("CTX 2\n"), TEMPLATE(T, ctx_print)(ctx2), flint_printf("\n");
                flint_printf("Embed\n"),
                    TEMPLATE(B, mat_print_pretty)(embed), flint_printf("\nProject\n"),
                    TEMPLATE(B, mat_print_pretty)(project), flint_printf("\nComposition\n"),
                    TEMPLATE(B, mat_print_pretty)(comp), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(B, mat_clear)(embed);
            TEMPLATE(B, mat_clear)(project);
            TEMPLATE(B, mat_clear)(comp);
            TEMPLATE(B, mat_clear)(one);
            TEMPLATE(B, poly_clear)(minpoly);
            TEMPLATE(B, poly_clear)(modulus2);
            TEMPLATE(T, clear)(gen1, ctx1);
            TEMPLATE(T, ctx_clear)(ctx1);
            TEMPLATE(T, clear)(gen2, ctx2);
            TEMPLATE(T, ctx_clear)(ctx2);
        }
    }

    TEST_FUNCTION_END(state);
}
#endif
#endif
