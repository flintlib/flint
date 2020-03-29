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
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

int main(void)
{
    slong rep, r, c1, c2, nreps = 100;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A, B, C;
    nmod_sparse_mat_t window1, window2;
    FLINT_TEST_INIT(state);

    flint_printf("concat_horizontal....");
    fflush(stdout);


    for (rep = 0; rep < nreps; rep++)
    {
        r = n_randint(state, 200);
        c1 = n_randint(state, 200);
        c2 = n_randint(state, 200);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c1, mod);
        nmod_sparse_mat_init(B, r, c2, mod);
        nmod_sparse_mat_init(C, r, c1+c2, mod);

        nmod_sparse_mat_randtest(A, state, 0, c1);
        nmod_sparse_mat_randtest(B, state, 0, c2);
        nmod_sparse_mat_randtest(C, state, 0, c1+c2);

        nmod_sparse_mat_concat_horizontal(C, A, B);        
        
        nmod_sparse_mat_window_init(window1, C, 0, 0, r, c1);
        nmod_sparse_mat_window_init(window2, C, 0, c1, r, c1+c2);

        if (!(nmod_sparse_mat_equal(window1, A) && nmod_sparse_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("B = \n");
            nmod_sparse_mat_print_pretty(B);
            flint_printf("A | B = \n");
            nmod_sparse_mat_print_pretty(C);
            flint_printf("window1 = \n");
            nmod_sparse_mat_print_pretty(window1);
            flint_printf("window2 = \n");
            nmod_sparse_mat_print_pretty(window2);
            flint_printf("FAIL: results not equal\n");
            abort();
        }

        nmod_sparse_mat_window_clear(window1);
        nmod_sparse_mat_window_clear(window2);

        nmod_sparse_mat_init(window1, r, c1, mod);
        nmod_sparse_mat_init(window2, r, c2, mod);
        nmod_sparse_mat_split_horizontal(window1, window2, C, c1);

       if (!(nmod_sparse_mat_equal(window1, A) && nmod_sparse_mat_equal(window2, B)))
        {
            flint_printf("A = \n");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("B = \n");
            nmod_sparse_mat_print_pretty(B);
            flint_printf("A | B = \n");
            nmod_sparse_mat_print_pretty(C);
            flint_printf("window1 = \n");
            nmod_sparse_mat_print_pretty(window1);
            flint_printf("window2 = \n");
            nmod_sparse_mat_print_pretty(window2);
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
