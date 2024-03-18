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

TEST_FUNCTION_START(n_factor_ecm, state)
{
    int i, j, k, result, fails;
    mp_limb_t prime1, prime2, prod, f, mod;

    fails = 0;

    for (i = 10; i < 64; i += 5)
    {
        for (j = i; j < 64 - i; j += 5)
        {
            for (k = 0; k < flint_test_multiplier(); k++)
            {
                prime1 = n_randprime(state, i, 1);
                prime2 = n_randprime(state, j, 1);
                prod = prime1 * prime2;

                result = n_factor_ecm(&f, (i + j) << 2, 1000, 50000, state, prod);

                if (result)
                {
                    mod = prod % f;
                    if ((mod != 0) || (f == prod) || (f == 1))
                        TEST_FUNCTION_FAIL(
                                "Wrong answer from stage %d\n"
                                "Number: %wu = %wu * %wu\n"
                                "Factor found: %wu\n",
                                result, prod, prime1, prime2, f);
                }
                else
                    fails += 1;
            }
        }
    }

    if (fails > 2 * flint_test_multiplier())
        TEST_FUNCTION_FAIL("Too many unsuccessful factorizations, %d\n", fails);

    TEST_FUNCTION_END(state);
}
