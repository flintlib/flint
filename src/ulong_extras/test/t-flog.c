/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_flog, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        ulong a = 0, b = 0, k, x;

        while (a < 1)
            a = n_randtest(state);
        while (b < 2)
            b = n_randtest(state);

        k = n_flog(a, b);
        x = n_pow(b, k);

        result = (x <= a && a / b < x);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "a = %wu\n"
                    "b = %wu\n"
                    "x = %wu\n"
                    "k = %wu\n",
                    a, b, x, k);
    }

    TEST_FUNCTION_END(state);
}
