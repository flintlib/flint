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
    TEMPLATE(T, sparse_vec_t) u, v;
    slong *inds;
    TEMPLATE(T, struct) *vals; 
    FLINT_TEST_INIT(state);
    
    flint_printf("construction from elements....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        TEMPLATE(T, sparse_vec_init) (u, ctx);
        TEMPLATE(T, sparse_vec_init) (v, ctx);
        TEMPLATE(T, sparse_vec_randtest) (u, state, nnz, len, ctx);
        TEMPLATE(T, sparse_vec_randtest) (v, state, nnz, len, ctx);
        
        /* Construct v from entries of u */
        inds = flint_malloc(nnz * sizeof(*inds));
        vals = flint_malloc(nnz * sizeof(*vals));
        for (i = 0; i < nnz; ++i) 
        {
            inds[i] = u->entries[i].ind;
            TEMPLATE(T, init) (&vals[i], ctx);
            TEMPLATE(T, set) (&vals[i], u->entries[i].val, ctx);
        }
        TEMPLATE(T, sparse_vec_from_entries) (v, inds, vals, nnz, ctx);

        if (!TEMPLATE(T, sparse_vec_equal) (u, v, 0, ctx))
        {
            flint_printf("FAIL: u != v\n");
            abort();
        }
        flint_free(inds);
        for (i = 0; i < nnz; ++i) 
        {
            TEMPLATE(T, clear) (&vals[i], ctx);
        }
        flint_free(vals);
        TEMPLATE(T, sparse_vec_clear) (u, ctx);
        TEMPLATE(T, sparse_vec_clear) (v, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif