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
    slong rep, len, nnz, i;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, t) *val;
    TEMPLATE(T, sparse_vec_t) u, v;
    TEMPLATE(T, struct) *w, *x;
    FLINT_TEST_INIT(state);
    
    flint_printf("conversion to/from dense vector....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        TEMPLATE(T, sparse_vec_init) (u, ctx);
        TEMPLATE(T, sparse_vec_init) (v, ctx);
        w = _TEMPLATE(T, vec_init) (len, ctx);
        x = _TEMPLATE(T, vec_init) (len, ctx);

        TEMPLATE(T, sparse_vec_randtest) (u, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_randtest) (v, state, nnz, len, ctx);
        
        TEMPLATE(T, sparse_vec_to_dense) (w, u, len, ctx);
        TEMPLATE(T, sparse_vec_from_dense) (v, w, len, ctx);

        for (i = 0; i < len; ++i)
        {
            val = TEMPLATE(T, sparse_vec_at) (u, i, ctx);
            if ((val == NULL && (!TEMPLATE(T, is_zero) (&w[i], ctx))) || (val != NULL && !(TEMPLATE(T, equal) (*val, &w[i], ctx))))
            {
                flint_printf("FAIL: u[%wd] != v[%wd]\n", i, i);
                abort();
            }
        }
        if (!TEMPLATE(T, sparse_vec_equal) (u, v, 0, ctx))
        {
            flint_printf("FAIL: u != v\n");
            abort();
        }

        _TEMPLATE(T, vec_randtest) (w, state, len, ctx);
        TEMPLATE(T, sparse_vec_from_dense) (u, w, len, ctx);
        TEMPLATE(T, sparse_vec_to_dense) (x, u, len, ctx);

        if (!_TEMPLATE(T, vec_equal) (w, x, len, ctx))
        {
            flint_printf("FAIL: w != x\n");
            abort();
        }
        TEMPLATE(T, sparse_vec_clear) (u, ctx);
        TEMPLATE(T, sparse_vec_clear) (v, ctx);
        _TEMPLATE(T, vec_clear) (w, len, ctx);
        _TEMPLATE(T, vec_clear) (x, len, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif