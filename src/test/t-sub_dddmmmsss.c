/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(sub_dddmmmsss, state)
{
    int i, j, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        mp_limb_t s[3], t[3], a[3], b[3];

        for (j = 0; j < 3; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        sub_dddmmmsss(s[2], s[1], s[0], a[2], a[1], a[0], b[2], b[1], b[0]);

        mpn_sub_n(t, a, b, 3);

        result = ((s[2] == t[2]) && (s[1] == t[1]) && (s[0] == t[0]));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a[2] = %wu, a[1] = %wu, a[0] = %wu\n", a[2], a[1], a[0]);
            flint_printf("b[2] = %wu, b[1] = %wu, b[0] = %wu\n", b[2], b[1], b[0]);
            flint_printf("s[2] = %wu, s[1] = %wu, s[0] = %wu\n", s[2], s[1], s[0]);
            flint_printf("t[2] = %wu, t[1] = %wu, t[0] = %wu\n", t[2], t[1], t[0]);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
