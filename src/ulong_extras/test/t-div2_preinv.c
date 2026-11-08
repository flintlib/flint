/*
    Copyright (C) 2009, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_div2_preinv, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, n, q1, q2;

        d = n_randtest_not_zero(state);
        n = n_randtest(state);
        if (n_randint(state, 2) == 0)
        {
            /* close to a multiple of d */
            ulong s = n_randint(state, 3);
            n = (n / d) * d;
            if (n <= UWORD_MAX - s)
                n += s;
        }

        dinv = n_preinvert_limb(d);

        q1 = n_div2_preinv(n, d, dinv);
        q2 = n / d;

        result = (q1 == q2);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu, d = %wu, dinv = %wu\n"
                    "q1 = %wu, q2 = %wu\n",
                    n, d, dinv, q1, q2);
    }

    TEST_FUNCTION_END(state);
}
