/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_is_oddprime_small, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        ulong d;
        int type;

        type = n_randint(state, 2);

        if (type == 0)
        {
            /* Test that primes pass the test */
            int bits = 2 + n_randint(state,
                    FLINT_BIT_COUNT(FLINT_ODDPRIME_SMALL_CUTOFF - 1) - 1);

            do
                d = n_randprime(state, bits, 1);
            while (d > FLINT_ODDPRIME_SMALL_CUTOFF);

            result = n_is_oddprime_small(d);
        }
        else
        {
            /* Test that not too many composites pass */
            do
                d = n_randint(state, FLINT_ODDPRIME_SMALL_CUTOFF) | 1;
            while (n_is_probabprime(d));

            result = !n_is_oddprime_small(d);
        }

        if (!result)
            TEST_FUNCTION_FAIL(
                    "type %d\n"
                    "d = %wu\n",
                    type, d);
    }

    TEST_FUNCTION_END(state);
}
