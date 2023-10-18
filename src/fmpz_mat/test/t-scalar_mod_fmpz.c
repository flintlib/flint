/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_scalar_mod_fmpz, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, Amod;
        fmpz_t mod;

        slong rows, cols;

        rows = n_randint(state, 20);
        cols = n_randint(state, 20);

        fmpz_mat_init(A, rows, cols);
        fmpz_mat_init(Amod, rows, cols);
        fmpz_init(mod);

        fmpz_randtest_not_zero(mod, state, 100);

        fmpz_mat_randtest(A, state, 100);
        fmpz_mat_randtest(Amod, state, 100);

        fmpz_mat_scalar_mod_fmpz(Amod, A, mod);
        fmpz_mat_scalar_mod_fmpz(A, A, mod);

        if (!fmpz_mat_equal(A, Amod))
        {
            flint_printf("FAIL: aliasing!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(Amod);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
