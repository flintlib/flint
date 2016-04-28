/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main(void)
{
    fmpz_t prime1, prime2, primeprod, fac, modval;
    int i, j, k, fails;

    FLINT_TEST_INIT(state);

    fmpz_init(prime1);
    fmpz_init(prime2);
    fmpz_init(primeprod);
    fmpz_init(fac);
    fmpz_init(modval);

    fails = 0;

    flint_printf("ecm....");
    fflush(stdout);

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
                    abort();
                }
            }
        }
    }

    if (fails > flint_test_multiplier())
    {
        printf("FAIL : ECM failed too many times (%d times)\n", fails);
        abort();
    }

    fmpz_clear(prime1);
    fmpz_clear(prime2);
    fmpz_clear(primeprod);
    fmpz_clear(fac);
    fmpz_clear(modval);
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
