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

int
main(void)
{
    slong rep, len, nnz, i;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_vec_t u, v;
    slong *inds;
    mp_limb_t *vals;
    FLINT_TEST_INIT(state);
    
    flint_printf("construction from elements....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);

        nmod_sparse_vec_init(u);
        nmod_sparse_vec_init(v);
        nmod_sparse_vec_randtest(u, state, nnz, len, mod);
        nmod_sparse_vec_randtest(v, state, nnz, len, mod);
        
        /* Construct v from entries of u */
        inds = flint_malloc(nnz * sizeof(*inds));
        vals = flint_malloc(nnz * sizeof(*vals));
        for (i = 0; i < nnz; ++i) 
        {
            inds[i] = u->entries[i].ind;
            vals[i] = u->entries[i].val;
        }
        nmod_sparse_vec_from_entries(v, inds, vals, nnz);

        if (!nmod_sparse_vec_equal(u, v, 0))
        {
            flint_printf("FAIL: u != v\n");
            abort();
        }
        flint_free(inds);
        flint_free(vals);
        nmod_sparse_vec_clear(u);
        nmod_sparse_vec_clear(v);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
