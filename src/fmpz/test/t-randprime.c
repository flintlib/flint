/*
    Authored 2015 by Daniel S. Roche; US Government work in the public domain.

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_randprime, state)
{
    int i;
    flint_bitcnt_t bits;

    for (bits=2; bits < 36; ++bits)
    {
        for (i=0; i < 1 * flint_test_multiplier(); ++i)
        {
            fmpz_t p;

            fmpz_init(p);

            fmpz_randprime(p, state, bits, 1);

            if (fmpz_bits(p) != bits || !_fmpz_is_canonical(p))
            {
                flint_printf("FAIL: not %wu bits\n", bits);
                fmpz_print(p); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (fmpz_is_prime(p) != 1)
            {
                flint_printf("FAIL: not prime, %wu bits\n", bits);
                fmpz_print(p); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(p);
        }
    }

    for (; bits < 130; ++bits)
    {
        /* at this point, the chances of a collision are less than
         * one in a billion. So it should never happen.
         */

        for (i=0; i < 1 + flint_test_multiplier()/5; ++i)
        {
            int j;
            fmpz p[2];

            fmpz_init(p+0);
            fmpz_init(p+1);

            for (j=0; j<2; ++j)
            {
                fmpz_randprime(p+j, state, bits, 0);

                if (fmpz_bits(p+j) != bits || !_fmpz_is_canonical(p+j))
                {
                    flint_printf("FAIL: not %wu bits\n", bits);
                    fmpz_print(p+j); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }

                if (fmpz_is_prime(p+j) == 0)
                {
                    flint_printf("FAIL: not prime, %wu bits\n", bits);
                    fmpz_print(p+j); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            if (fmpz_equal(p+0, p+1))
            {
                flint_printf("FAIL: returned a duplicate\n");
                fmpz_print(p+0); flint_printf("\n");
                fmpz_print(p+1); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(p+0);
            fmpz_clear(p+1);
        }
    }

    TEST_FUNCTION_END(state);
}
