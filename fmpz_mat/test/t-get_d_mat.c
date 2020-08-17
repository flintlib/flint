/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
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
#include "fmpz.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);


    flint_printf("get_d_mat....");
    fflush(stdout);

    /* set entries of an fmpz_mat, convert to d_mat and then check that
       the entries remain same */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        d_mat_t B;
        slong j, k;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        d_mat_init(B, rows, cols);

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                fmpz_set_ui(fmpz_mat_entry(A, j, k), 3 * j + 7 * k);
            }
        }

        fmpz_mat_get_d_mat(B, A);

        for (j = 0; j < rows; j++)
        {
            for (k = 0; k < cols; k++)
            {
                if (d_mat_entry(B, j, k) != 3 * j + 7 * k)
                {
                    flint_printf("FAIL: j = %wd, k = %wd\n", j, k);
                    fmpz_mat_print_pretty(A);
                    d_mat_print(B);
                    abort();
                }
            }
        }

        fmpz_mat_clear(A);
        d_mat_clear(B);
    }



    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
