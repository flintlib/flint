/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("window_init/clear....");
    fflush(stdout);


    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
    	TEMPLATE(T, ctx_t) ctx;

    	TEMPLATE(T, mat_t) a, w;
        slong j, r1, r2, c1, c2;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, mat_init) (a, rows, cols, ctx);
        TEMPLATE(T, mat_randtest) (a, state, ctx);

        r2 = n_randint(state, rows + 1);
        c2 = n_randint(state, cols + 1);
        if (r2)
            r1 = n_randint(state, r2);
        else
            r1 = 0;
        if (c2)
            c1 = n_randint(state, c2);
        else
            c1 = 0;

        TEMPLATE(T, mat_window_init) (w, a, r1, c1, r2, c2, ctx);

        for (j = 0; j < r2 - r1; j++)
        	 _TEMPLATE(T, vec_zero) (w->rows[j], c2 - c1, ctx);

        TEMPLATE(T, mat_window_clear) (w, ctx);
        TEMPLATE(T, mat_clear) (a, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }


    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
