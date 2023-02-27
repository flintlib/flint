/*
    Copyright (C) 2010 Fredrik Johansson

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
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    fmpz_mat_t A;
    slong i, m, n, b, d, r;

    FLINT_TEST_INIT(state);
 
    flint_printf("rank....");
    fflush(stdout);

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 35);
        n = n_randint(state, 35);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);
            fmpz_mat_init(A, m, n);
            fmpz_mat_randrank(A, state, r, b);
            if (r != fmpz_mat_rank(A))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }
            fmpz_mat_clear(A);
        }
    }

    /* Dense */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 35);
        n = n_randint(state, 35);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);
            fmpz_mat_init(A, m, n);
            fmpz_mat_randrank(A, state, r, b);
            fmpz_mat_randops(A, state, d);
            if (r != fmpz_mat_rank(A))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }
            fmpz_mat_clear(A);
        }
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
