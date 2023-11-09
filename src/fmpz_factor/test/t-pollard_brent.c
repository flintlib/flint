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
#include "fmpz.h"
#include "fmpz_factor.h"

TEST_FUNCTION_START(fmpz_factor_pollard_brent, state)
{
    fmpz_t prime1, prime2, prime3, prime4, primeprod, fac, modval;
    int i, j, k, fails;

    fmpz_init(prime1);
    fmpz_init(prime2);
    fmpz_init(prime3);
    fmpz_init(prime4);
    fmpz_init(primeprod);
    fmpz_init(fac);
    fmpz_init(modval);

    fails = 0;

    for (i = 10; i < 26; i += 5)
    {
        for (j = 0; j < 10 * flint_test_multiplier(); j++)
        {
            fmpz_set_ui(prime1, n_randprime(state, i, 1));
            fmpz_set_ui(prime2, n_randprime(state, i, 1));
            fmpz_set_ui(prime3, n_randprime(state, i, 1));
            fmpz_set_ui(prime4, n_randprime(state, i, 1));

            fmpz_mul(prime1, prime1, prime2);
            fmpz_mul(prime3, prime3, prime4);
            fmpz_mul(primeprod, prime1, prime3);

            k = fmpz_factor_pollard_brent(fac, state, primeprod, 5, 2500);

            if (k == 0)
                fails += 1;
            else
            {
                fmpz_mod(modval, primeprod, fac);
                k = fmpz_cmp_ui(modval, 0);
                if (k != 0)
                {
                    printf("FAIL : Wrong factor calculated\n");
                    printf("n : ");
                    fmpz_print(primeprod);
                    printf(" factor calculated : ");
                    fmpz_print(fac);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
    }

    if (fails > 2*flint_test_multiplier())
    {
        printf("FAIL : Pollard Rho failed too many times (%d times)\n", fails);
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(prime1);
    fmpz_clear(prime2);
    fmpz_clear(prime3);
    fmpz_clear(prime4);
    fmpz_clear(primeprod);
    fmpz_clear(fac);
    fmpz_clear(modval);

    TEST_FUNCTION_END(state);
}
