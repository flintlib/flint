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
    

    flint_printf("scalar_add/submul_nmod_mat_ui....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        nmod_mat_t M;
        slong rows, cols;
        ulong c;
        ulong mod;

        rows = n_randint(state, 10);
        cols = n_randint(state, 10);
        mod = n_randtest_prime(state, 0);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(B, rows, cols);
        fmpz_mat_init(C, rows, cols);
        nmod_mat_init(M, rows, cols, mod);
        c = n_randtest(state);

        nmod_mat_randtest(M, state);
        fmpz_mat_set_nmod_mat_unsigned(A, M);

        fmpz_mat_randtest(B, state, 100);
        fmpz_mat_set(C, B);

        fmpz_mat_scalar_addmul_nmod_mat_ui(B, M, c);
        fmpz_mat_scalar_addmul_ui(C, A, c);

        if (!fmpz_mat_equal(B, C))
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        nmod_mat_clear(M);
    }

    

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
