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
    slong rep, len, nnz;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_vec_t u, v, w, x;

    FLINT_TEST_INIT(state);
    
    flint_printf("add/sub....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        len = n_randint(state, 40);
        nnz = n_randint(state, len+1);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);

        nmod_sparse_vec_init(u);
        nmod_sparse_vec_init(v);
        nmod_sparse_vec_init(w);
        nmod_sparse_vec_init(x);

        nmod_sparse_vec_randtest(u, state, nnz, len, mod);
        nmod_sparse_vec_randtest(v, state, nnz, len, mod);

        nmod_sparse_vec_add(w, u, v, mod);
        nmod_sparse_vec_sub(x, w, v, mod);

        if (!nmod_sparse_vec_equal(u, x, 0))
        {
            flint_printf("FAIL: u != u+v-v\n");
            abort();
        }

        nmod_sparse_vec_add(u, u, v, mod);
        if (!nmod_sparse_vec_equal(u, w, 0))
        {
            flint_printf("FAIL: (u += v) != u + v\n");
            abort();
        }

        nmod_sparse_vec_sub(u, u, v, mod);
        if (!nmod_sparse_vec_equal(u, x, 0))
        {
            flint_printf("FAIL: ((u += v) -= v) != u+v-v\n");
            abort();
        }

        nmod_sparse_vec_add(u, v, u, mod);
        if (!nmod_sparse_vec_equal(u, w, 0))
        {
            flint_printf("FAIL: (u = v + u) != u+v\n");
            abort();
        }

        nmod_sparse_vec_sub(v, u, v, mod);
        if (!nmod_sparse_vec_equal(v, x, 0))
        {
            flint_printf("FAIL: (u = v + u) != u+v\n");
            abort();
        }
        nmod_sparse_vec_clear(u);
        nmod_sparse_vec_clear(v);
        nmod_sparse_vec_clear(w);
        nmod_sparse_vec_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
