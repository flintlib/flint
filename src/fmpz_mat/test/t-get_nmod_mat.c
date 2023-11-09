/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_get_nmod_mat, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A;
        nmod_mat_t M, M2;
        slong rows, cols;
        mp_limb_t mod;

        rows = n_randint(state, 50);
        cols = n_randint(state, 50);

        mod = n_randtest_prime(state, 0);

        nmod_mat_init(M, rows, cols, mod);
        nmod_mat_init(M2, rows, cols, mod);
        fmpz_mat_init(A, rows, cols);

        nmod_mat_randtest(M, state);

        if (i % 2 == 0)
            fmpz_mat_set_nmod_mat(A, M);
        else
            fmpz_mat_set_nmod_mat_unsigned(A, M);

        fmpz_mat_scalar_mul_ui(A, A, UWORD(2));
        nmod_mat_add(M, M, M);
        fmpz_mat_get_nmod_mat(M2, A);

        if (!nmod_mat_equal(M, M2))
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        nmod_mat_clear(M);
        nmod_mat_clear(M2);
    }

    TEST_FUNCTION_END(state);
}
