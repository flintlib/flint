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

TEST_FUNCTION_START(fmpz_factor_ecm, state)
{
    fmpz_t prime1, prime2, primeprod, fac, modval;
    int i, j, k, fails;

    fmpz_init(prime1);
    fmpz_init(prime2);
    fmpz_init(primeprod);
    fmpz_init(fac);
    fmpz_init(modval);

    fails = 0;

    for (i = 35; i <= 50; i += 5)
    {
        for (j = 0; j < flint_test_multiplier(); j++)
        {
            fmpz_set_ui(prime1, n_randprime(state, i, 1));
            fmpz_set_ui(prime2, n_randprime(state, i, 1));

            fmpz_mul(primeprod, prime1, prime2);

            k = fmpz_factor_ecm(fac, i << 2, 2000, 50000, state, primeprod);

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

    if (fails > flint_test_multiplier())
    {
        printf("FAIL : ECM failed too many times (%d times)\n", fails);
        fflush(stdout);
        flint_abort();
    }

    /* Tests for hangs and crashes, don't care about result */

#if FLINT64
    fmpz_set_ui(prime1, 1123047674690129);
    fmpz_set_ui(prime2, 66049336315331);
    fmpz_mul(primeprod, prime1, prime2);
    fmpz_mul(primeprod, primeprod, primeprod);
    fmpz_factor_ecm(fac, 1, 100, 1000, state, primeprod);

    /* (p*q)^2 for p and q of 53 bits */
    for (i = 0; i < 5; i++)
    {
        fmpz_set_ui(prime1, n_randprime(state, 53, 1));
	fmpz_set_ui(prime2, n_randprime(state, 53, 1));
	fmpz_mul(primeprod, prime1, prime2);
	fmpz_mul(primeprod, primeprod, primeprod);
	fmpz_factor_ecm(fac, 212, 2000, 50000, state, primeprod);
    }

    /* p^2*q*r for p and q of 53 bits */
    for (i = 0; i < 5; i++)
    {
        fmpz_set_ui(primeprod, n_randprime(state, 53, 1));
        fmpz_mul(primeprod, primeprod, primeprod);
        fmpz_mul_ui(primeprod, primeprod, n_randprime(state, 53, 1));
        fmpz_mul_ui(primeprod, primeprod, n_randprime(state, 53, 1));
        fmpz_factor_ecm(fac, 212, 2000, 50000, state, primeprod);
    }
#endif

    fmpz_clear(prime1);
    fmpz_clear(prime2);
    fmpz_clear(primeprod);
    fmpz_clear(fac);
    fmpz_clear(modval);

    TEST_FUNCTION_END(state);
}
