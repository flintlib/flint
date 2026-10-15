/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_randprime_mod, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong bits, m, r, p;

        bits = 2 + n_randint(state, FLINT_BITS - 1);
        m = 1 + n_randint(state, FLINT_MIN(UWORD(1) << (bits - 1), 1000));
        do
            r = n_randint(state, m);
        while (n_gcd(r, m) != 1);

        /* a short range may contain no prime in the progression */
        if (bits < 16)
        {
            ulong n, lo = UWORD(1) << (bits - 1), hi = (UWORD(1) << bits) - 1;
            int exists = 0;
            for (n = lo; n <= hi && !exists; n++)
                exists = (n % m == r && n_is_prime(n));
            if (!exists)
                continue;
        }

        p = n_randprime_mod(state, bits, r, m, n_randint(state, 2));

        if (!n_is_prime(p) || p % m != r || FLINT_BIT_COUNT(p) != bits)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits = %wu, r = %wu, m = %wu, p = %wu\n", bits, r, m, p);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
