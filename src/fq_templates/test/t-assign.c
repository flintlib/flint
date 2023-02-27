/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);

    flint_printf("assignments... ");
    fflush(stdout);

    /* Check that gen does not segfault */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) x;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        
        TEMPLATE(T, init)(x, ctx);

        TEMPLATE(T, gen)(x, ctx);

        TEMPLATE(T, clear)(x, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
