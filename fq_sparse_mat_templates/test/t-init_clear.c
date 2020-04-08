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
    slong rep, r, c, i;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);

        if (!TEMPLATE(T, sparse_mat_is_zero) (A, ctx))
        {
            flint_printf("FAIL: A not zero!\n");
            abort();
        }
        for (i = 0; i < r; i++)
        {
            if (!TEMPLATE(T, sparse_vec_is_zero) (&A->rows[i], ctx))
            {
                flint_printf("FAIL: row %wd not zero!\n", i);
                abort();
            }
        }

        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
