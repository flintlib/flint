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

TEST_FUNCTION_START(n_mod2_preinv, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        ulong d, dinv, n, r1, r2, r3;

        d = n_randtest_not_zero(state);
        if (n_randint(state, 8) == 0)
            d = UWORD(1) << n_randint(state, FLINT_BITS);

        n = n_randtest(state);
        if (n_randint(state, 2) == 0)
        {
            /* close to a multiple of d, near the top of the word */
            ulong q = UWORD_MAX / d - n_randint(state, 3);
            ulong s = n_randint(state, 4);
            n = q * d;
            if (s == 3)
                s = d - 1;
            if (n <= UWORD_MAX - s)
                n += s;
        }

        dinv = n_preinvert_limb(d);

        r1 = n_mod2_preinv(n, d, dinv);
        r2 = n % d;
        r3 = n_mod2_preinv_lazy(n, d, dinv);

        /* the lazy result is r2 or r2 + d, and below d + 2 */
        result = (r1 == r2) && (r3 == r2 || (r2 <= 1 && r3 == r2 + d));
        if (!result)
            TEST_FUNCTION_FAIL(
                    "n = %wu, d = %wu, dinv = %wu\n"
                    "r1 = %wu, r2 = %wu, r3 = %wu\n",
                    n, d, dinv, r1, r2, r3);
    }

    TEST_FUNCTION_END(state);
}
