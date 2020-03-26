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
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, mod, i, j, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        nmod_sparse_mat_t A;

        m = n_randint(state, 50);
        do
            mod = n_randtest_not_zero(state);
        while(mod <= 1);

        nmod_sparse_mat_init(A, m, mod);
        if(A->nnz != UWORD(0)) {
            flint_printf("FAIL: nnz not zero!\n");
            abort();
        }
        if(A->entries != NULL) {
            flint_printf("FAIL: entries not null!\n");
            abort();
        }
        if(A->c != UWORD(0)) {
            flint_printf("FAIL: c not 0!\n");
            abort();
        }
        
        for (i = 0; i < m; i++)
        {
            if (A->row_starts[i] != UWORD(0))
            {
                flint_printf("FAIL: row start not zero!\n");
                abort();
            }
            if (A->row_nnz[i] != UWORD(0))
            {
                flint_printf("FAIL: row nnz not zero!\n");
                abort();
            }
        }

        if (A->mod.n != mod)
        {
            flint_printf("FAIL: bad modulus\n");
            abort();
        }

        nmod_sparse_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
