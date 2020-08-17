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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("scalar_mul/divexact_fmpz....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t C, B, A;
        slong rows, cols;
        fmpz_t c;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);
        fmpz_init(c);

        fmpz_randtest(c, state, 100);
        fmpz_mat_randtest(A, state, 100);

        fmpz_mat_scalar_mul_fmpz(B, A, c);

        if (fmpz_is_zero(c))
        {
            if (!fmpz_mat_is_zero(B))
            {
                flint_printf("FAIL!\n");
                abort();
            }
        }
        else
        {
            fmpz_mat_scalar_divexact_fmpz(C, B, c);

            if (!fmpz_mat_equal(C, A))
            {
                flint_printf("FAIL!\n");
                abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_clear(c);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
