/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

int
main(void)
{
    int i;
    
    FLINT_TEST_INIT(state);

    flint_printf("embed matrices... ");
    fflush(stdout);

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
            abort();
        }

        TEMPLATE(B, mat_clear)(embed);
        TEMPLATE(B, mat_clear)(project);
        TEMPLATE(B, mat_clear)(one);
        TEMPLATE(T, clear)(gen, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check random isomorphism */
    /* Check proper embeddings */
    /* Check Artin-Schreier case */

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


#endif
#endif
