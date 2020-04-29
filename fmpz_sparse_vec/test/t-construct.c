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

int
main(void)
{
    slong rep, bits, len, nnz, i;
    fmpz_sparse_vec_t u, v;
    slong *inds;
    fmpz *vals;
    FLINT_TEST_INIT(state);
    
    flint_printf("construction from elements....");
    fflush(stdout);

    for (rep = 0; rep < 1; rep++)
    {
        do bits = n_randint(state, 100);
        while (bits < UWORD(2));
        len = n_randint(state, 10);
        nnz = n_randint(state, len+1);

        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);
        
        /* Construct v from entries of u */
        inds = flint_malloc(nnz * sizeof(*inds));
        vals = flint_malloc(nnz * sizeof(*vals));
        for (i = 0; i < nnz; ++i) 
        {
            fmpz_init(&vals[i]);
            inds[i] = u->entries[i].ind;
            fmpz_set(&vals[i], u->entries[i].val);
        }
        fmpz_sparse_vec_from_entries(v, inds, vals, nnz);

        if (!fmpz_sparse_vec_equal(u, v, 0))
        {
            fmpz_sparse_vec_print_pretty(u, 0, 0);
            fmpz_sparse_vec_print_pretty(v, 0, 0);
            abort();
        }
        flint_free(inds);
        for (i = 0; i < nnz; ++i)
        {
            fmpz_clear(&vals[i]);
        }
        flint_free(vals);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
