/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_init_clear, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t a;
        mp_limb_t mod;
        slong j, k;
        slong rows = n_randint(state, 100);
        slong cols = n_randint(state, 100);
        mod = n_randtest_prime(state, 0);

        nmod_poly_mat_init(a, rows, cols, mod);

        for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
                nmod_poly_zero(nmod_poly_mat_entry(a, j, k));

        nmod_poly_mat_clear(a);
    }

    TEST_FUNCTION_END(state);
}
