/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_window_init_clear, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t a, w;
        slong j, k, r1, r2, c1, c2;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        fmpz_mat_init(a, rows, cols);
        fmpz_mat_randtest(a, state, n_randint(state, 200) + 1);

        r2 = n_randint(state, rows + 1);
        c2 = n_randint(state, cols + 1);
        if (r2)
            r1 = n_randint(state, r2);
        else
            r1 = 0;
        if (c2)
            c1 = n_randint(state, c2);
        else
            c1 = 0;

        fmpz_mat_window_init(w, a, r1, c1, r2, c2);

        for (j = 0; j < r2 - r1; j++)
            for (k = 0; k < c2 - c1; k++)
                fmpz_zero(fmpz_mat_entry(w, j, k));

        fmpz_mat_window_clear(w);
        fmpz_mat_clear(a);
    }

    TEST_FUNCTION_END(state);
}
