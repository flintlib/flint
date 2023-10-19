/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_init_clear, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t a;
        slong j, k;
        slong rows = n_randint(state, 100);
        slong cols = n_randint(state, 100);
        fmpz_t mod;

        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 200);
        fmpz_abs(mod, mod);

        fmpz_mod_mat_init(a, rows, cols, mod);

        for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
                fmpz_zero(fmpz_mod_mat_entry(a, j, k));

        fmpz_mod_mat_clear(a);
        fmpz_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
