/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_lehman, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        ulong n1, n2;
        int type;

        type = n_randint(state, 10);

        if (type == 0)
        {
            /* Test random products of two primes */
            int bits1, bits2;
            bits1 = 2 + n_randint(state, FLINT_MIN(FLINT_BITS, 53) - 3);
            bits2 = 2 + n_randint(state, FLINT_MIN(FLINT_BITS, 53) - 1 - bits1);

            n1 = n_randprime(state, bits1, 1) * n_randprime(state, bits2, 1);
        }
        else
        {
            /* Test random numbers */
            do
            {
                int bits = 2 + n_randint(state, FLINT_MIN(FLINT_BITS, 53) - 1);
                n1 = n_randtest_bits(state, bits);
            } while (n_is_prime(n1));
        }

#if FLINT64
        /* Test a specific bug (FIXME: What bug?) */
        if (ix == 0)
            n1 = UWORD(72528697) * UWORD(73339073);
#endif

        n2 = n_factor_lehman(n1);

        result = ((n1 % n2) == UWORD(0) && n1 != n2);
        if (!result)
            TEST_FUNCTION_FAIL("n1 = %wu, n2 = %wu\n", n1, n2);
    }

    TEST_FUNCTION_END(state);
}
