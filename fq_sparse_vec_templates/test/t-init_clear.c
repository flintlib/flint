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

static void check_zero(TEMPLATE(T, sparse_vec_t) vec)
{
    if (vec->nnz != UWORD(0))
{
        flint_printf("FAIL: nnz not zero!\n");
        abort();
    }
    if (vec->entries != NULL)
{
        flint_printf("FAIL: entries not null!\n");
        abort();
    }
}

int
main(void)
{
    slong rep, len, nnz, i;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_vec_t) vec;
    FLINT_TEST_INIT(state);
    
    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, sparse_vec_init) (vec, ctx);

        check_zero(vec);

        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);

        TEMPLATE(T, sparse_vec_randtest) (vec, state, nnz, len, ctx);
        
        if (nnz == 0) check_zero(vec);
        else 
        {
            for (i = 0; i < nnz; ++i)
            {
                TEMPLATE(T, sparse_entry_struct) *e = &vec->entries[i];
                if (e->ind >= len)
                {
                    flint_printf("FAIL: found index %wd >= %wd!\n", e->ind, len);
                    abort();
                }
                if (TEMPLATE(T, is_zero) (e->val, ctx))
                {
                    flint_printf("FAIL: found 0 value\n");
                    abort();
                }
                if (i > 0 && e->ind <= e[-1].ind)
                {
                    flint_printf("FAIL: found index %wd <= previous index %wd\n", e->ind, e[-1].ind);
                    abort();
                }
            }
        }

        TEMPLATE(T, sparse_vec_clear) (vec, ctx);
        check_zero(vec);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif