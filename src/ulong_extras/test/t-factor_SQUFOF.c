/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_SQUFOF, state)
{
    int i, result;
    ulong count = UWORD(0);

    for (i = 0; i < 300 * flint_test_multiplier(); i++) /* Test random numbers */
    {
        mp_limb_t n1, n2;

        do
        {
            n1 = n_randtest_bits(state, n_randint(state, FLINT_BITS) + 1);
        } while (n_is_prime(n1) || (n1 < UWORD(2)));

#if FLINT64
        n2 = n_factor_SQUFOF(n1, 10000);
#else
        n2 = n_factor_SQUFOF(n1, 2000);
#endif

        if (n2)
        {
            count++;
            result = ((n1%n2) == UWORD(0));
            if (!result)
                TEST_FUNCTION_FAIL("n1 = %wu, n2 = %wu\n", n1, n2);
        }
    }

    if (count < 280 * flint_test_multiplier())
        TEST_FUNCTION_FAIL("Only %wu numbers factored\n", count);

    TEST_FUNCTION_END(state);
}
