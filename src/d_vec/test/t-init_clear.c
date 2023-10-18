/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "d_vec.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(d_vec_init_clear, state)
{
    int i;

    /* check if memory management works properly */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        double *a;
        slong j, len = n_randint(state, 100) + 1;

        a = _d_vec_init(len);
        for (j = 0; j < len; j++)
            a[j] = 0;

        _d_vec_clear(a);
    }

    TEST_FUNCTION_END(state);
}
