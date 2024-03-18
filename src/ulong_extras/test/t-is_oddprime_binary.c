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

TEST_FUNCTION_START(n_is_oddprime_binary, state)
{
    slong ix;
    slong cutoff = 100000;
    int result;

    for (ix = 0; ix < 20000 * flint_test_multiplier(); ix++)
    {
        ulong d;
        int type;

        type = n_randint(state, 2);

        /* NOTE: n_is_oddprime_binary(n) requires n to be odd and n > 16 */
        if (type == 0)
        {
            /* Test that primes pass the test */
            int bits = 5 + n_randint(state, 13);
            do
                d = n_randprime(state, bits, 1);
            while (d > cutoff);

            result = n_is_oddprime_binary(d);
        }
        else
        {
            /* Test that not too many composites pass */
            do
                d = (17 + n_randint(state, cutoff - 17)) | 1;
            while (n_is_probabprime(d));

            result = !n_is_oddprime_binary(d);
        }

        if (!result)
            TEST_FUNCTION_FAIL(
                    "type %d\n"
                    "d = %wu\n",
                    type, d);
    }

    TEST_FUNCTION_END(state);
}
