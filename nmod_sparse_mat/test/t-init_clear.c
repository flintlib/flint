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
    slong rep, r, c, i;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A;
    FLINT_TEST_INIT(state);
    

    flint_printf("init/clear....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        r = n_randint(state, 200);
        c = n_randint(state, 200);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);

        if (!nmod_sparse_mat_is_zero(A))
        {
            flint_printf("FAIL: A not zero!\n");
            abort();
        }
        for (i = 0; i < r; i++)
        {
            if (!nmod_sparse_vec_is_zero(&A->rows[i])) 
            {
                flint_printf("FAIL: row %wd not zero!\n", i);
                abort();
            }
        }

        nmod_sparse_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
