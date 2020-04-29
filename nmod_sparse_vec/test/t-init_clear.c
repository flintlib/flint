/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "ulong_extras.h"

static void check_zero(nmod_sparse_vec_t vec)
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
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_vec_t vec;
    FLINT_TEST_INIT(state);
    
    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);

        nmod_sparse_vec_init(vec);
        check_zero(vec);

        nmod_sparse_vec_randtest(vec, state, nnz, len, mod);
        
        if (nnz == 0) check_zero(vec);
        else 
        {
            for (i = 0; i < nnz; ++i)
            {
                nmod_sparse_entry_struct *e = &vec->entries[i];
                if (e->ind >= len)
                {
                    flint_printf("FAIL: found index %wd >= %wd!\n", e->ind, len);
                    abort();
                }
                if (e->val == UWORD(0) || e->val >= n)
                {
                    flint_printf("FAIL: found value %wd (not in (0,%wd))\n", e->val, n);
                    abort();
                }
                if (i > 0 && e->ind <= e[-1].ind)
                {
                    flint_printf("FAIL: found index %wd <= previous index %wd\n", e->ind, e[-1].ind);
                    abort();
                }
            }
        }

        nmod_sparse_vec_clear(vec);
        check_zero(vec);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
