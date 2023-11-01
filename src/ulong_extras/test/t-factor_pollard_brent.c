/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_factor_pollard_brent, state)
{
    mp_limb_t prime1, prime2, primeprod, fac, modval;
    int i, j, k, l, fails;

    fails = 0;

    for (l = 5; l < 26; l += 5)
    {
        for (i = l; i < 26 && i + l <= FLINT_BITS; i += 5)
        {
            for (j = 0; j < 10 * flint_test_multiplier(); j++)
            {
                do {
                    prime1 = n_randtest_bits(state, l);
                    prime2 = n_randtest_bits(state, i);
                    primeprod = prime1 * prime2;
                } while (primeprod < 1);

                k = n_factor_pollard_brent(&fac, state, primeprod, 5, 2500);

                if (k == 0)
                    fails += 1;
                else
                {
                    modval = primeprod % fac;
                    if (modval != 0)
                    {
                        flint_printf("FAIL : Wrong factor calculated\n");
                        flint_printf("n : %wu\n", primeprod);
                        flint_printf("Factor calculated: %wn\n", fac);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }
    }

#if FLINT64
    if (fails > flint_test_multiplier())
#else
    if (fails > 2 * flint_test_multiplier())
#endif
    {
        printf("FAIL : Pollard Rho - Brent failed too many times (%d times)\n", fails);
        fflush(stdout);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
