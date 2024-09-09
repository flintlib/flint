/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "test_helpers.h"

TEST_FUNCTION_START(sub_dddmmmsss, state)
{
    int i, j, n, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong s[8], t[8], a[8], b[8];
        int aliasing;

        for (j = 0; j < 8; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        aliasing = n_randint(state, 2);

        for (n = 2; n < 8; n++)
        {
            if (aliasing)
            {
                for (j = 0; j < 8; j++)
                    s[j] = a[j];

                if (n == 2)
                    NN_SUB_2(s, s, b);
                else if (n == 3)
                    NN_SUB_3(s, s, b);
                else if (n == 4)
                    NN_SUB_4(s, s, b);
                else if (n == 5)
                    NN_SUB_5(s, s, b);
                else if (n == 6)
                    NN_SUB_6(s, s, b);
                else if (n == 7)
                    NN_SUB_7(s, s, b);
                else if (n == 8)
                    NN_SUB_8(s, s, b);
            }
            else
            {
                if (n == 2)
                    NN_SUB_2(s, a, b);
                else if (n == 3)
                    NN_SUB_3(s, a, b);
                else if (n == 4)
                    NN_SUB_4(s, a, b);
                else if (n == 5)
                    NN_SUB_5(s, a, b);
                else if (n == 6)
                    NN_SUB_6(s, a, b);
                else if (n == 7)
                    NN_SUB_7(s, a, b);
                else if (n == 8)
                    NN_SUB_8(s, a, b);
            }

            mpn_sub_n(t, a, b, n);

            result = flint_mpn_equal_p(s, t, n);

            if (!result)
            {
                TEST_FUNCTION_FAIL(
                        "Aliasing: %d\n"
                        "n = %d\n"
                        "a = %{ulong*}\n"
                        "b = %{ulong*}\n"
                        "s = %{ulong*}\n"
                        "t = %{ulong*}\n",
                        aliasing,
                        n,
                        a, n,
                        b, n,
                        s, n,
                        t, n);
            }
        }
    }

    TEST_FUNCTION_END(state);
}
