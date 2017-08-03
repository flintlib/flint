/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, e, d;
    FLINT_TEST_INIT(state);

    flint_printf("embed... ");
    fflush(stdout);

    /* Check isomorphism to self */
    for (i = 0; i < 4 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, t) a, b;

        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, init) (a, ctx);
        TEMPLATE(T, init) (b, ctx);

        TEMPLATE(T, embed_gens) (a, b, ctx, ctx);
        
        d = TEMPLATE(T, ctx_degree)(ctx);
        for (e = 0; e < d; e++) {
            if (TEMPLATE(T, equal)(a, b, ctx))
                break;
            TEMPLATE(T, frobenius)(b, b, 1, ctx);
        }
        
        if (e == d)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), TEMPLATE(T, ctx_print)(ctx), flint_printf("\n");
            flint_printf("a: "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b: "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear) (a, ctx);
        TEMPLATE(T, clear) (b, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
