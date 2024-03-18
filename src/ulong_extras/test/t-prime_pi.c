/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_prime_pi, state)
{
    int n;

    for (n=1; n<10000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        if ((n_prime_pi(n-1)+1 == n_prime_pi(n)) != n_is_prime(n))
            TEST_FUNCTION_FAIL("expected pi(%d) + 1 = pi(%d)\n", n-1, n);
    }

    for (n=1; n<5000 * FLINT_MIN(10, flint_test_multiplier()); n++)
    {
        if (n_prime_pi(n_nth_prime(n)) != n)
            TEST_FUNCTION_FAIL("expected pi(prime(%d)) = %d\n", n, n);
    }

    TEST_FUNCTION_END(state);
}
