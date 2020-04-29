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
#include "fmpz_sparse_vec.h"
#include "ulong_extras.h"

int main(void)
{
    slong rep, bits, len, nnz;
    fmpz_sparse_vec_t u, v, w;
    fmpz_sparse_vec_t window1, window2;
    FLINT_TEST_INIT(state);

    flint_printf("concat....");
    fflush(stdout);

    for (rep = 0; rep < 100; rep++)
    {
        do bits = n_randint(state, 200);
        while (bits < UWORD(2));
        len = n_randint(state, 200);
        nnz = n_randint(state, len+1);

        fmpz_sparse_vec_init(u);
        fmpz_sparse_vec_init(v);
        fmpz_sparse_vec_init(w);

        fmpz_sparse_vec_randtest(u, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(v, state, nnz, len, bits);
        fmpz_sparse_vec_randtest(w, state, nnz, len, bits);

         fmpz_sparse_vec_concat(w, u, v, len);        
        
        fmpz_sparse_vec_window_init(window1, w, 0, len);
        fmpz_sparse_vec_window_init(window2, w, len, 2*len);

        if (!(fmpz_sparse_vec_equal(window1, u, 0) && fmpz_sparse_vec_equal(window2, v, len)))
        {
            flint_printf("u = ");
            fmpz_sparse_vec_print_pretty(u, 0, len);
            flint_printf("v = \n");
            fmpz_sparse_vec_print_pretty(v, 0, len);
            flint_printf("u | v = \n");
            fmpz_sparse_vec_print_pretty(w, 0, len);
            flint_printf("window1 = \n");
            fmpz_sparse_vec_print_pretty(window1, 0, len);
            flint_printf("window2 = \n");
            fmpz_sparse_vec_print_pretty(window2, len, len);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        fmpz_sparse_vec_window_clear(window1);
        fmpz_sparse_vec_window_clear(window2);
        
        fmpz_sparse_vec_init(window1);
        fmpz_sparse_vec_init(window2);
        fmpz_sparse_vec_split(window1, window2, w, len);
        if (!(fmpz_sparse_vec_equal(window1, u, 0) && fmpz_sparse_vec_equal(window2, v, 0)))
        {
            flint_printf("u = ");
            fmpz_sparse_vec_print_pretty(u, 0, len);
            flint_printf("v = \n");
            fmpz_sparse_vec_print_pretty(v, 0, len);
            flint_printf("u | v = \n");
            fmpz_sparse_vec_print_pretty(w, 0, len);
            flint_printf("window1 = \n");
            fmpz_sparse_vec_print_pretty(window1, 0, len);
            flint_printf("window2 = \n");
            fmpz_sparse_vec_print_pretty(window2, 0, len);
            flint_printf("FAIL: results not equal\n");
            abort();
        }
        fmpz_sparse_vec_clear(window1);
        fmpz_sparse_vec_clear(window2);
        fmpz_sparse_vec_clear(u);
        fmpz_sparse_vec_clear(v);
        fmpz_sparse_vec_clear(w);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
