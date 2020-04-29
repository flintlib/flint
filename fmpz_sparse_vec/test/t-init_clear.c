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
#include "fmpz_sparse_vec.h"
#include "ulong_extras.h"

static void check_zero(fmpz_sparse_vec_t vec)
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
    slong rep, len, nnz, bits, i;
    fmpz_sparse_vec_t vec;
    FLINT_TEST_INIT(state);
    
    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        len = n_randint(state, 50);
        nnz = n_randint(state, len+1);
        fmpz_sparse_vec_init(vec);
        check_zero(vec);

        fmpz_sparse_vec_randtest(vec, state, nnz, len, bits);
        
        if (nnz == 0) check_zero(vec);
        else 
        {
            for (i = 0; i < nnz; ++i)
            {
                fmpz_sparse_entry_struct *e = &vec->entries[i];
                if (e->ind >= len)
                {
                    flint_printf("FAIL: found index %wd >= %wd!\n", e->ind, len);
                    abort();
                }
                if (fmpz_is_zero(e->val))
                {
                    flint_printf("FAIL: found zero value\n");
                    abort();
                }
                if (i > 0 && e->ind <= e[-1].ind)
                {
                    flint_printf("FAIL: found index %wd <= previous index %wd\n", e->ind, e[-1].ind);
                    abort();
                }
            }
        }

        fmpz_sparse_vec_clear(vec);
        check_zero(vec);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
