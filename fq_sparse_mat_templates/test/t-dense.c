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
    slong rep, r, c;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, B;
    TEMPLATE(T, mat_t) C, D;
    FLINT_TEST_INIT(state);
    

    flint_printf("conversion to/from dense matrix....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r = n_randint(state, 10);
        c = n_randint(state, 10);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (B, r, c, ctx);
        TEMPLATE(T, mat_init) (C, r, c, ctx);
        TEMPLATE(T, mat_init) (D, r, c, ctx);
        
        TEMPLATE(T, sparse_mat_randtest) (A, state, 0, c, ctx);
        TEMPLATE(T, sparse_mat_to_dense) (C, A, ctx);
        TEMPLATE(T, sparse_mat_from_dense) (B, C, ctx);
        
        if (!TEMPLATE(T, sparse_mat_equal) (A, B, ctx))
        {
            flint_printf("FAIL: A != B\n");
            abort();
        }

        TEMPLATE(T, mat_randtest) (C, state, ctx);
        TEMPLATE(T, sparse_mat_from_dense) (A, C, ctx);
        TEMPLATE(T, sparse_mat_to_dense) (D, A, ctx);
        
        if (!TEMPLATE(T, mat_equal) (C, D, ctx))
        {
            flint_printf("FAIL: C != D\n");
            abort();
        }
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (B, ctx);
        TEMPLATE(T, mat_clear) (C, ctx);
        TEMPLATE(T, mat_clear) (D, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
