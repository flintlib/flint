/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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
        slong i;
        mp_limb_t p, q;

        n_primes_init(iter);
        q = 0;
        for (i = 0; i < 200000; i++)
        {
            p = n_primes_next(iter);
            q = n_nextprime(q, 0);

            if (p != q)
            {
                flint_printf("FAIL\n");
                flint_printf("i = %wu, p = %wu, q = %wu\n", i, p, q);
                fflush(stdout);
                flint_abort();
            }
        }

        n_primes_clear(iter);
    }

    /* count primes */
    for (n = 0; n < 10; n++)
    {
        n_primes_t iter;
        mp_limb_t s, p, r;

        const unsigned int primepi[10] = {
            0, 4, 25, 168, 1229, 9592, 78498, 664579, 5761455, 50847534
        };

        r = n_pow(10, n);

        n_primes_init(iter);
        s = 0;
        while ((p = n_primes_next(iter)) <= r)
            s++;

        if (s != primepi[n])
        {
            flint_printf("FAIL\n");
            flint_printf("pi(10^%wd) = %u, computed = %wu\n", n, primepi[n], s);
            fflush(stdout);
            flint_abort();
        }

        n_primes_clear(iter);
    }

    TEST_FUNCTION_END(state);
}
