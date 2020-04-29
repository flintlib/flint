/*
    Copyright (C) 2010 Fredrik Johansson

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

int
main(void)
{
    slong rep, r, c, nreps = 1000;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, B, C;
    FLINT_TEST_INIT(state);
    

    flint_printf("transpose....");
    fflush(stdout);

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < nreps; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r = n_randint(state, 20);
        c = n_randint(state, 20);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (B, c, r, ctx);
        TEMPLATE(T, sparse_mat_init) (C, r, c, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, 0, c, ctx);
        TEMPLATE(T, sparse_mat_randtest) (B, state, 0, r, ctx);
        TEMPLATE(T, sparse_mat_randtest) (C, state, 0, c, ctx);

        TEMPLATE(T, sparse_mat_transpose) (B, A, ctx);
        TEMPLATE(T, sparse_mat_transpose) (C, B, ctx);
        
        if (!TEMPLATE(T, sparse_mat_equal) (C, A, ctx))
        {
            flint_printf("FAIL: C != A\n");
            abort();
        }

        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (B, ctx);
        TEMPLATE(T, sparse_mat_clear) (C, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
