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

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("embed... ");
    fflush(stdout);

    /* Check isomorphism to self */
    for (i = 0; i < 4 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, ev_a, ev_b;
        TEMPLATE(B, poly_t) minpoly;
        TEMPLATE(T, poly_t) minpoly_fq;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, init) (a, ctx);
        TEMPLATE(T, init) (b, ctx);
        TEMPLATE(T, init) (ev_a, ctx);
        TEMPLATE(T, init) (ev_b, ctx);
        TEMPLATE(B, poly_init) (minpoly, 
                                TEMPLATE(B, poly_modulus)(TEMPLATE(T, ctx_modulus)(ctx)));
        TEMPLATE(T, poly_init) (minpoly_fq, ctx);

        TEMPLATE(T, embed_gens) (a, b, minpoly, ctx, ctx);

        TEMPLATE4(T, poly_set, B, poly)(minpoly_fq, minpoly, ctx);
        TEMPLATE3(T, poly_evaluate, T)(ev_a, minpoly_fq, a, ctx);
        TEMPLATE3(T, poly_evaluate, T)(ev_b, minpoly_fq, b, ctx);
        
        if (!TEMPLATE(T, is_zero)(ev_a, ctx) || !TEMPLATE(T, is_zero)(ev_b, ctx))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("a: "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b: "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("minpoly: "), TEMPLATE(B, poly_print_pretty)(minpoly, "x"), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear) (a, ctx);
        TEMPLATE(T, clear) (b, ctx);
        TEMPLATE(T, clear) (ev_a, ctx);
        TEMPLATE(T, clear) (ev_b, ctx);
        TEMPLATE(B, poly_clear) (minpoly);
        TEMPLATE(T, poly_clear) (minpoly_fq, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


#endif
#endif
