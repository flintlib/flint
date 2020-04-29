/*
    Copyright (C) 2011 Fredrik Johansson

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
    slong rep, r, c;
    TEMPLATE(T, t) a, cc;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, B, C, D;
    FLINT_TEST_INIT(state);
    
    flint_printf("scalar_mul....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, init) (a, ctx);
        TEMPLATE(T, init) (cc, ctx);
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        TEMPLATE(T, randtest) (cc, state, ctx);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (B, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (C, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (D, r, c, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, 0, c, ctx);

        TEMPLATE(T, TEMPLATE(sparse_mat_scalar_mul, T)) (B, A, a, ctx);
        TEMPLATE(T, one) (cc, ctx);
        TEMPLATE(T, sub) (cc, a, cc, ctx);
        TEMPLATE(T, TEMPLATE(sparse_mat_scalar_mul, T)) (C, A, cc, ctx);

        /* c*A - (c-1)*A == A */
        TEMPLATE(T, sparse_mat_sub) (D, B, C, ctx);

        if (!TEMPLATE(T, sparse_mat_equal) (A, D, ctx))
        {
            flint_printf("FAIL\n");
            abort();
        }

         /* Aliasing */
        TEMPLATE(T, TEMPLATE(sparse_mat_scalar_mul, T)) (C, A, a, ctx);
        TEMPLATE(T, TEMPLATE(sparse_mat_scalar_mul, T)) (A, A, a, ctx);

        if (!TEMPLATE(T, sparse_mat_equal) (A, C, ctx))
        {
            flint_printf("FAIL\n");
            abort();
        }

        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (B, ctx);
        TEMPLATE(T, sparse_mat_clear) (C, ctx);
        TEMPLATE(T, sparse_mat_clear) (D, ctx);
        TEMPLATE(T, clear) (a, ctx);
        TEMPLATE(T, clear) (cc, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
