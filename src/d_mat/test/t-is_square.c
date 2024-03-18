/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "d_mat.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(d_mat_is_square, state)
{
    int i;

    /* check if a non-zero value is returned for non-square matrices */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        d_mat_t A;
        slong rows = n_randint(state, 10);
        slong cols = n_randint(state, 10);

        d_mat_init(A, rows, cols);

        if (d_mat_is_square(A) != (rows == cols))
            TEST_FUNCTION_FAIL("");

        d_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
