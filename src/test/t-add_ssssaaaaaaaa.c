/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(add_ssssaaaaaaaa, state)
{
    int i, j, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        mp_limb_t s[4], t[4], a[4], b[4];
        int aliasing;

        for (j = 0; j < 4; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        aliasing = n_randint(state, 2);

        if (aliasing)
        {
            for (j = 0; j < 4; j++)
                s[j] = a[j];

            add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], s[3], s[2], s[1], s[0],
                                                     b[3], b[2], b[1], b[0]);
        }
        else
        {
            add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], a[3], a[2], a[1], a[0],
                                                     b[3], b[2], b[1], b[0]);
        }

        mpn_add_n(t, a, b, 4);

        result = ((s[3] == t[3]) && (s[2] == t[2]) && (s[1] == t[1]) && (s[0] == t[0]));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "Aliasing: %d\n"
                    "a[3] = %wu, a[2] = %wu, a[1] = %wu, a[0] = %wu\n"
                    "b[3] = %wu, b[2] = %wu, b[1] = %wu, b[0] = %wu\n"
                    "s[3] = %wu, s[2] = %wu, s[1] = %wu, s[0] = %wu\n"
                    "t[3] = %wu, t[2] = %wu, t[1] = %wu, t[0] = %wu\n",
                    aliasing,
                    a[3], a[2], a[1], a[0],
                    b[3], b[2], b[1], b[0],
                    s[3], s[2], s[1], s[0],
                    t[3], t[2], t[1], t[0]);
    }

    TEST_FUNCTION_END(state);
}
