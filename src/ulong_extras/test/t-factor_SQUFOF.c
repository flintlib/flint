/*
    Copyright (C) 2009 William Hart
                  2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

#define ITER_SMALL     10
#define ITER_BIG     1000
#define BITS_SMALL    2 + RND(FLINT_BITS / 2 - 1)
#define BITS_BIG      FLINT_BITS / 2 + RND(FLINT_BITS / 2) + 1

TEST_FUNCTION_START(n_factor_SQUFOF, state)
{
    int i, result;
    slong nr_tests = 1000 * FLINT_MAX(1, flint_test_multiplier()),
          nr_success = 0;

    for (i = 0; i < nr_tests; i++)
    {
        ulong num, fac;

        if (n_randint(state, 4))
        {
            do num = RNDBITS(BITS_SMALL);
            while (n_is_prime(num));
            fac = n_factor_SQUFOF(num, ITER_SMALL);
        }
        else
        {
            do num = RNDBITS(BITS_BIG);
            while (n_is_prime(num));
            fac = n_factor_SQUFOF(num, ITER_BIG);
        }

        if (fac)
        {
            nr_success++;
            result = (num % fac == 0);

            if (!result)
                TEST_FUNCTION_FAIL("num = %wu, fac = %wu\n", num, fac);
        }
    }

    if (nr_success < 0.95 * nr_tests)
        TEST_FUNCTION_FAIL("Only %wu / %wd = %f worked\n",
                           nr_success, nr_tests,
                           (double) nr_success / nr_tests);

    TEST_FUNCTION_END(state);
}

#undef ITER_SMALL
#undef ITER_BIG
#undef BITS_SMALL
#undef BITS_BIG
