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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_sparse_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong rep, bits, r, c, nreps = 1000, num_bad_fail = 0;
    fmpz_sparse_mat_t A;
    fmpz_mat_t X, B, AX, AXm, Bm;
    fmpz_t mod, det;
    int success;

    FLINT_TEST_INIT(state);

    flint_printf("solve_dixon....");
    fflush(stdout);    

    for (rep = 0; rep < nreps; rep++)
    {
        do bits = n_randint(state, 100);
        while (bits < UWORD(2));
        r = n_randint(state, 20);
        c = n_randint(state, 20);

        fmpz_sparse_mat_init(A, r, r);
        fmpz_mat_init(B, r, c);
        fmpz_mat_init(X, r, c);
        fmpz_init(mod);

        fmpz_sparse_mat_randtest(A, state, 0, r, bits);
        fmpz_mat_randtest(B, state, 1+n_randint(state, 2)*bits);
        
        success = fmpz_sparse_mat_solve_dixon(X, mod, A, B);
        if (!success)
        {
            fmpz_init(det);
            fmpz_sparse_mat_det(det, A);
            if (!fmpz_is_zero(det)) 
            {
                num_bad_fail++;
            }
            fmpz_clear(det);
        }
        else
        {
            fmpz_mat_init(AX, r, c);
            fmpz_mat_init(AXm, r, c);
            fmpz_mat_init(Bm, r, c);

            fmpz_sparse_mat_mul_mat(AX, A, X);
            fmpz_mat_scalar_mod_fmpz(AXm, AX, mod);
            fmpz_mat_scalar_mod_fmpz(Bm, B, mod);

            if (!fmpz_mat_equal(AXm, Bm) || !success)
            {
                flint_printf("FAIL:\n");
                flint_printf("AX != B!\n");
                flint_printf("A:\n"),      fmpz_sparse_mat_print_pretty(A),  flint_printf("\n");
                flint_printf("B:\n"),      fmpz_mat_print_pretty(B),  flint_printf("\n");
                flint_printf("X:\n"),      fmpz_mat_print_pretty(X),  flint_printf("\n");
                flint_printf("mod = "),    fmpz_print(mod),           flint_printf("\n");
                flint_printf("AX:\n"),     fmpz_mat_print_pretty(AX), flint_printf("\n");
                abort();
            }
            fmpz_mat_clear(Bm);
            fmpz_mat_clear(AX);
            fmpz_mat_clear(AXm);
        }
        

        fmpz_sparse_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    if (num_bad_fail != 0) flint_printf("Number of unexpected failures: %wd\n", num_bad_fail);
    return 0;
}
