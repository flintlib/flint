/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_primes, state)
{
    slong n;

    /* compare with n_nextprime */
    {
        n_primes_t iter;
        mp_limb_t p, q;

        n_primes_init(iter);
        q = 0;
        for (n = 0; n < 100000 * flint_test_multiplier(); n++)
        {
            p = n_primes_next(iter);
            q = n_nextprime(q, 0);

            if (p != q)
                TEST_FUNCTION_FAIL("n = %wu, p = %wu, q = %wu\n", n, p, q);
        }

        n_primes_clear(iter);
    }

    /* count primes */
    {
        ulong r, s, p;
        n_primes_t iter;
        slong n_max = flint_test_multiplier() > 1.0 ? 10 : flint_test_multiplier() < 1.0 ? 8 : 9;
        const unsigned int primepi[10] =
        {
            0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534
        };

        n_primes_init(iter);

        for (n = 0, r = 1, s = 0; n < n_max; n++, r *= 10, s++)
        {
            while ((p = n_primes_next(iter)) <= r)
                s++;

            if (s != primepi[n])
                TEST_FUNCTION_FAIL("pi(10^%wd) = %u, computed = %wu\n", n, primepi[n], s);
        }

        n_primes_clear(iter);
    }

    TEST_FUNCTION_END(state);
}
