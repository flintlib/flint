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
    mp_ptr w, x;
    FLINT_TEST_INIT(state);
    
    flint_printf("conversion to/from dense vector....");
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
        w = _nmod_vec_init(len);
        x = _nmod_vec_init(len);

        nmod_sparse_vec_randtest(u, state, nnz, len, mod);
        nmod_sparse_vec_randtest(v, state, nnz, len, mod);
        
        nmod_sparse_vec_to_dense(w, u, len);
        nmod_sparse_vec_from_dense(v, w, len);

        for (i = 0; i < len; ++i)
        {
            mp_limb_t *val = nmod_sparse_vec_at(u, i);
            if ((w[i] == 0 && val != NULL) || (w[i] != 0 && (val == NULL || *val != w[i])))
        {
                flint_printf("FAIL: u[%wd] != v[%wd]\n", i, i);
                abort();
            }
        }
        if (!nmod_sparse_vec_equal(u, v, 0))
        {
            flint_printf("FAIL: u != v\n");
            abort();
        }

        _nmod_vec_randtest(w, state, len, mod);
        nmod_sparse_vec_from_dense(u, w, len);
        nmod_sparse_vec_to_dense(x, u, len);

        if (!_nmod_vec_equal(w, x, len))
        {
            flint_printf("FAIL: w != x\n");
            abort();
        }
        nmod_sparse_vec_clear(u);
        nmod_sparse_vec_clear(v);
        _nmod_vec_clear(w);
        _nmod_vec_clear(x);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
