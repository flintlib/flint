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
#include "nmod_vec.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int main(void)
{
    slong rep, r1, r2, c, nreps = 100;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A, B, C;
    nmod_sparse_mat_t window1, window2;
    FLINT_TEST_INIT(state);


    flint_printf("concat_vertical....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        r1 = n_randint(state, 100);
        r2 = n_randint(state, 100);
        c = n_randint(state, 100);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r1, c, mod);
        nmod_sparse_mat_init(B, r2, c, mod);
        nmod_sparse_mat_init(C, r1+r2, c, mod);

        nmod_sparse_mat_randtest(A, state, 0, c);
        nmod_sparse_mat_randtest(B, state, 0, c);
        nmod_sparse_mat_randtest(C, state, 0, c);

        nmod_sparse_mat_concat_vertical(C, A, B);
        
        nmod_sparse_mat_window_init(window1, C, 0, 0, r1, c);
        nmod_sparse_mat_window_init(window2, C, r1, 0, r1+r2, c);

        if (!(nmod_sparse_mat_equal(window1, A) && nmod_sparse_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("B = \n");
            nmod_sparse_mat_print_pretty(B);
            flint_printf("A concat_vertical B = \n");
            nmod_sparse_mat_print_pretty(C);
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        nmod_sparse_mat_window_clear(window1);
        nmod_sparse_mat_window_clear(window2);

        nmod_sparse_mat_init(window1, r1, c, mod);
        nmod_sparse_mat_init(window2, r2, c, mod);
        nmod_sparse_mat_split_vertical(window1, window2, C, r1);
        
        if (!(nmod_sparse_mat_equal(window1, A) && nmod_sparse_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("B = \n");
            nmod_sparse_mat_print_pretty(B);
            flint_printf("A concat_vertical B = \n");
            nmod_sparse_mat_print_pretty(C);
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        nmod_sparse_mat_window_clear(window1);
        nmod_sparse_mat_window_clear(window2);
        
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_sparse_mat_clear(C);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
