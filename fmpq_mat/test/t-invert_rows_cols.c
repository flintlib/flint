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
#include "fmpq.h"
#include "fmpq_mat.h"
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
        fmpq_mat_t A, B;

        n = n_randint(state, 10);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);

        fmpq_mat_randtest(A, state, 1+n_randint(state, 3));

        fmpq_mat_set(B, A);
        fmpq_mat_invert_rows(A, NULL);
        fmpq_mat_invert_cols(A, NULL);

        for (i = 0; i < A->r; i++)
        {
            for (j =0; j < A->c; j++)
            {
                if (fmpq_cmp(fmpq_mat_entry(B, i, j), fmpq_mat_entry(A, A->r - i - 1, A->c - j - 1)) != 0)
                {
                    flint_printf("FAIL: B != A\n");
                    flint_printf("A:\n");
                    fmpq_mat_print(A);
                    flint_printf("B:\n");
                    fmpq_mat_print(B);
                    abort();
                }
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
