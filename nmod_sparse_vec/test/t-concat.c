/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "ulong_extras.h"

int main(void)
{
    slong rep, len, nnz;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_vec_t u, v, w;
    nmod_sparse_vec_t window1, window2;
    FLINT_TEST_INIT(state);

    flint_printf("concat....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);

        nmod_sparse_vec_init(u);
        nmod_sparse_vec_init(v);
        nmod_sparse_vec_init(w);

        nmod_sparse_vec_randtest(u, state, nnz, len, mod);
        nmod_sparse_vec_randtest(v, state, nnz, len, mod);
        nmod_sparse_vec_randtest(w, state, nnz, len, mod);

        nmod_sparse_vec_concat(w, u, v, len);        
        
        nmod_sparse_vec_window_init(window1, w, 0, len);
        nmod_sparse_vec_window_init(window2, w, len, 2*len);

        if (!(nmod_sparse_vec_equal(window1, u, 0) && nmod_sparse_vec_equal(window2, v, len)))
        {
            flint_printf("u = ");
            nmod_sparse_vec_print_pretty(u, 0, len, mod);
            flint_printf("v = \n");
            nmod_sparse_vec_print_pretty(v, 0, len, mod);
            flint_printf("u | v = \n");
            nmod_sparse_vec_print_pretty(w, 0, len, mod);
            flint_printf("window1 = \n");
            nmod_sparse_vec_print_pretty(window1, 0, len, mod);
            flint_printf("window2 = \n");
            nmod_sparse_vec_print_pretty(window2, len, len, mod);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        nmod_sparse_vec_window_clear(window1);
        nmod_sparse_vec_window_clear(window2);

        nmod_sparse_vec_init(window1);
        nmod_sparse_vec_init(window2);
        nmod_sparse_vec_split(window1, window2, w, len);
        if (!(nmod_sparse_vec_equal(window1, u, 0) && nmod_sparse_vec_equal(window2, v, 0)))
        {
            flint_printf("u = ");
            nmod_sparse_vec_print_pretty(u, 0, len, mod);
            flint_printf("v = \n");
            nmod_sparse_vec_print_pretty(v, 0, len, mod);
            flint_printf("u | v = \n");
            nmod_sparse_vec_print_pretty(w, 0, len, mod);
            flint_printf("window1 = \n");
            nmod_sparse_vec_print_pretty(window1, 0, len, mod);
            flint_printf("window2 = \n");
            nmod_sparse_vec_print_pretty(window2, 0, len, mod);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        nmod_sparse_vec_window_clear(window1);
        nmod_sparse_vec_window_clear(window2);
        nmod_sparse_vec_clear(u);
        nmod_sparse_vec_clear(v);
        nmod_sparse_vec_clear(w);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
