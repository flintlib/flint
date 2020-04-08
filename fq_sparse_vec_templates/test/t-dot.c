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
    slong rep, len, nnz;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, t) a, b;
    TEMPLATE(T, sparse_vec_t) u, v;
    TEMPLATE(T, struct) *w, *x;
    FLINT_TEST_INIT(state);
    
    flint_printf("dot product....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, sparse_vec_init) (u, ctx);
        TEMPLATE(T, sparse_vec_init) (v, ctx);
        TEMPLATE(T, init) (a, ctx);
        TEMPLATE(T, init) (b, ctx);

        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        w = _TEMPLATE(T, vec_init) (len, ctx);
        x = _TEMPLATE(T, vec_init) (len, ctx);

        TEMPLATE(T, sparse_vec_randtest) (u, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_randtest) (v, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_to_dense) (w, u, len, ctx);
        TEMPLATE(T, sparse_vec_to_dense) (x, v, len, ctx);

        TEMPLATE(T, sparse_vec_dot) (a, u, v, ctx);
        _TEMPLATE(T, vec_dot) (b, w, x, len, ctx);

        if (!TEMPLATE(T, equal) (a, b, ctx))
        {
            flint_printf("Fail: sparse dot sparse\n");
            abort();
        }

        TEMPLATE(T, sparse_vec_dot_dense) (a, u, x, ctx);

        if (!TEMPLATE(T, equal) (a, b, ctx))
        {
            flint_printf("Fail: sparse dot dense\n");
            abort();
        }

        TEMPLATE(T, sparse_vec_clear) (u, ctx);
        TEMPLATE(T, sparse_vec_clear) (v, ctx);
        _TEMPLATE(T, vec_clear) (w, len, ctx);
        _TEMPLATE(T, vec_clear) (x, len, ctx);
        TEMPLATE(T, clear) (a, ctx);
        TEMPLATE(T, clear) (b, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif