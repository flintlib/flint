/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(add_sssaaaaaa, state)
{
    int i, j, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong s[3], t[3], a[3], b[3];

        for (j = 0; j < 3; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        add_sssaaaaaa(s[2], s[1], s[0], a[2], a[1], a[0], b[2], b[1], b[0]);

        mpn_add_n(t, a, b, 3);

        result = ((s[2] == t[2]) && (s[1] == t[1]) && (s[0] == t[0]));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "a[2] = %wu, a[1] = %wu, a[0] = %wu\n"
                    "b[2] = %wu, b[1] = %wu, b[0] = %wu\n"
                    "s[2] = %wu, s[1] = %wu, s[0] = %wu\n"
                    "t[2] = %wu, t[1] = %wu, t[0] = %wu\n",
                    a[2], a[1], a[0],
                    b[2], b[1], b[0],
                    s[2], s[1], s[0],
                    t[2], t[1], t[0]);
    }

    TEST_FUNCTION_END(state);
}
