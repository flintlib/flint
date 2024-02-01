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

TEST_FUNCTION_START(n_sqrtmod, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random integers */
    {
        mp_limb_t a, b, p, pinv;

        p = n_randtest_prime(state, 0);
        a = n_randtest(state) % p;

        b = n_sqrtmod(a, p);
        pinv = n_preinvert_limb(p);

        result = (b == 0 || n_mulmod2_preinv(b, b, p, pinv) == a);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "p = %wu\n"
                    "a = %wu\n"
                    "b = %wu\n",
                    p, a, b);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++) /* Test random squares */
    {
        mp_limb_t a, b, p, pinv;

        p = n_randtest_prime(state, 0);

        do
            b = n_randtest(state) % p;
        while (b == 0);

        pinv = n_preinvert_limb(p);
        a = n_mulmod2_preinv(b, b, p, pinv);

        b = n_sqrtmod(a, p);

        result = (n_mulmod2_preinv(b, b, p, pinv) == a);
        if (!result)
            TEST_FUNCTION_FAIL(
                    "p = %wu\n"
                    "a = %wu\n"
                    "b = %wu\n",
                    p, a, b);
    }

    TEST_FUNCTION_END(state);
}
