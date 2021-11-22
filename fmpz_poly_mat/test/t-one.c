/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

int
main(void)
{
    int iter;

    FLINT_TEST_INIT(state);

    flint_printf("one/is_one....");
    fflush(stdout);    

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_mat_t A;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_randtest(A, state, n_randint(state, 5),
            n_randint(state, 100));
        fmpz_poly_mat_one(A);

        if (!fmpz_poly_mat_is_one(A))
        {
            flint_printf("FAIL: expected matrix to be one\n");
            fflush(stdout);
            flint_abort();
        }

        if (m > 0 && n > 0)
        {
            m = n_randint(state, m);
            n = n_randint(state, n);

            if (m != n)
                fmpz_poly_randtest_not_zero(fmpz_poly_mat_entry(A, m, n),
                    state, 5, 5);
            else
                do { fmpz_poly_randtest_not_zero(fmpz_poly_mat_entry(A, m, n),
                    state, 5, 5); }
                while (fmpz_poly_is_one(fmpz_poly_mat_entry(A, m, n)));

            if (fmpz_poly_mat_is_one(A))
            {
                flint_printf("FAIL: expected matrix not to be one\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_poly_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
