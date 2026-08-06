/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_binvert, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong a, b, c;

        a = n_randtest(state) | 1;
        b = n_binvert(a);
        c = a * b;

        if (c != 1)
            TEST_FUNCTION_FAIL("a = %wu, b = %wu, c = %wu\n", a, b, c);
    }

    TEST_FUNCTION_END(state);
}
