/*
    Copyright (C) 2019 Tommy Hofmann Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong n, rep;
    FLINT_TEST_INIT(state);

    flint_printf("invert_rows/cols....");
    fflush(stdout);

    /* Rectangular inversion of rows and cols */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        slong i, j;
        fmpz_mat_t A, B;

        n = n_randint(state, 10);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);

        fmpz_mat_randtest(A, state, 1+n_randint(state, 3));

        fmpz_mat_set(B, A);
        fmpz_mat_invert_rows(A, NULL);
        fmpz_mat_invert_cols(A, NULL);

        for (i = 0; i < A->r; i++)
        {
            for (j =0; j < A->c; j++)
            {
                if (fmpz_cmp(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, A->r - i - 1, A->c - j - 1)) != 0)
                {
                    flint_printf("FAIL: B != A\n");
                    flint_printf("A:\n");
                    fmpz_mat_print_pretty(A);
                    flint_printf("B:\n");
                    fmpz_mat_print_pretty(B);
                    abort();
                }
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
