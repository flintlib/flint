/*
    Copyright (C) 2012 William Hart
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

/* Composite strong pseudoprimes from
   https://oeis.org/A014233 and https://doi.org/10.1090/mcom/3134 */
static const char * composites_is_prime[] = {
    "2047",
    "1373653",
    "25326001",
    "3215031751",
    "2152302898747",
    "3474749660383",
    "341550071728321",
    "230245660726188031",
    "3825123056546413051",
    "7395010240794120709381",
    "164280218643672633986221",
    "318665857834031151167461",
    "360681321802296925566181",
    "552727880697763694556181",
    "667636712015520329618581",
    "2995741773170734841812261",
    "3317044064679887385961981",
    "3404730287403079539471001",
    NULL,
};

TEST_FUNCTION_START(fmpz_is_prime, state)
{
    int i, result, r1;

    /* test table of composites */
    for (i = 0; composites_is_prime[i] != NULL; i++)
    {
        fmpz_t n;
        fmpz_init(n);
        fmpz_set_str(n, composites_is_prime[i], 10);

        r1 = fmpz_is_prime(n);
        if (r1 != 0)
        {
            flint_printf("FAIL:\n");
            fmpz_print(n); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
    }

    /* test primes always pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, F;

        fmpz_init(p);
        fmpz_init(F);

        do {
            fmpz_randbits(p, state, n_randint(state, 160) + 2);
            fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p));

        r1 = fmpz_is_prime(p);

        result = (r1 == 1 || r1 == -1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_clear(F);
    }

    /* test composites never pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, F;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(F);

        do {
            fmpz_randbits(p, state, n_randint(state, 80) + 2);
        } while (fmpz_cmp_ui(p, 2) < 0);
        do {
            fmpz_randbits(a, state, n_randint(state, 80) + 2);
        } while (fmpz_cmp_ui(a, 2) < 0);

        fmpz_mul(p, p, a);

        r1 = fmpz_is_prime(p);

        result = (r1 == 0 || r1 == -1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(F);
    }

    /* test issue 345 */
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_ui(t, 13567);
        fmpz_pow_ui(t, t, 145);

        if (fmpz_is_prime(t) != 0)
        {
            flint_printf("FAIL:\n");
            fmpz_print(t); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(t);
    }

    /* regression test */
    {
        fmpz_t p;
        fmpz_init_set_ui(p, 13567);
        fmpz_pow_ui(p, p, 145);
        result = fmpz_is_prime(p);
        if (result)
        {
            printf("FAIL\n");
            fmpz_print(p);
            printf("\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
