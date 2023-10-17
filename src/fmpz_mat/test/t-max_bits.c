/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_max_bits, state)
{
    slong m, n, rep, res1, res2;

    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_mat_init(A, m, n);
        fmpz_mat_randtest(A, state, 1 + n_randint(state, 100));

        res1 = fmpz_mat_max_bits(A);
        res2 = _fmpz_vec_max_bits(A->entries, m*n);

        if (res1 != res2)
        {
            flint_printf("FAIL!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
