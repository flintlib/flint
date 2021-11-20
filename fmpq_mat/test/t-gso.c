/*
    Copyright (C) 2014 Abhinav Baid

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
#include "fmpq_vec.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);


    flint_printf("gso....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C;
        fmpq_t dot;
        int j, k;

        slong m, n, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(C, n, m);

        fmpq_mat_randtest(A, state, bits);

        fmpq_mat_transpose(B, A);
        fmpq_mat_rref(B, B);

        fmpq_mat_gso(A, A);

        fmpq_mat_transpose(C, A);

        fmpq_init(dot);

        for (j = 0; j < n; j++)
        {
            for (k = j + 1; k < n; k++)
            {

                if (m)
                {
                   _fmpq_vec_dot(dot, C->rows[j], C->rows[k], m);

                   if (!fmpq_is_zero(dot))
                   {
                       flint_printf("FAIL:\n");
                       flint_printf("A:\n");
                       fmpq_mat_print(A);
                       fflush(stdout);
                       flint_abort();
                   }
                }
            }
        }

        fmpq_mat_rref(C, C);

        result = fmpq_mat_equal(B, C);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("B:\n");
            fmpq_mat_print(B);
            flint_printf("C:\n");
            fmpq_mat_print(C);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_clear(dot);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
