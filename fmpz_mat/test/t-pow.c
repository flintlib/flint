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

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("pow....");
    fflush(stdout);

    

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        slong i, n;
        ulong e;

        n = n_randint(state, 10);
        e = n_randint(state, 20);

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);
        fmpz_mat_init(C, n, n);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        fmpz_mat_pow(B, A, e);

        fmpz_mat_one(C);
        for (i = 0; i < e; i++)
            fmpz_mat_mul(C, C, A);

        if (!fmpz_mat_equal(C, B))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_pow(A, A, e);

        if (!fmpz_mat_equal(A, B))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
