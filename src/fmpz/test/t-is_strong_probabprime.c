/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

/* Composite strong pseudoprimes from https://oeis.org/A014233 */
static const char * composites_is_strong_probabprime[] = {
    "2047",
    "1373653",
    "25326001",
    "3215031751",
    "2152302898747",
    "3474749660383",
    "341550071728321",
    "341550071728321",
    "3825123056546413051",
    "3825123056546413051",
    "3825123056546413051",
    "318665857834031151167461",
    "3317044064679887385961981",
    NULL,
};

TEST_FUNCTION_START(fmpz_is_strong_probabprime, state)
{
    int i, result, count = 0;

    /* test table */
    {
        for (i = 0; composites_is_strong_probabprime[i] != NULL; i++)
        {
            int j;
            fmpz_t n, a;
            fmpz_init(n);
            fmpz_init(a);

            fmpz_set_str(n, composites_is_strong_probabprime[i], 10);

            for (j = 0; j <= i; j++)
            {
                fmpz_set_ui(a, n_nth_prime(j + 1));

                if (!fmpz_is_strong_probabprime(n, a))
                {
                    flint_printf("FAIL (composite expected to pass test):\n");
                    fmpz_print(n); printf("\n");
                    fmpz_print(a); printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_set_ui(a, n_nth_prime(i + 2
                    + (i == 6) + 2*(i==8) + (i==9))  /* because of repeated entries */
                );

            if (fmpz_is_strong_probabprime(n, a))
            {
                flint_printf("FAIL (composite expected to fail test):\n");
                fmpz_print(n); printf("\n");
                fmpz_print(a); printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(n);
            fmpz_clear(a);
        }
    }

    /* test primes always pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, b, F;

        fmpz_init(p);
        fmpz_init(b);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 330) + 2);
           fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p) || fmpz_cmp_ui(p, 2) == 0);

        do {
           fmpz_randbits(b, state, n_randint(state, 100) + 1);
           fmpz_abs(b, b);
        } while (fmpz_is_zero(b) || fmpz_is_one(b));

        result = fmpz_is_strong_probabprime(p, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            fmpz_print(b); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
        fmpz_clear(b);
        fmpz_clear(F);
    }

    /* test composites rarely pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, b, F;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(p, 2) < 0);
        do {
           fmpz_randbits(a, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(a, 2) < 0);

        do {
           fmpz_randbits(b, state, n_randint(state, 100) + 1);
           fmpz_abs(b, b);
        } while (fmpz_cmp_ui(b, 2) < 0);

        fmpz_mul(p, p, a);

        if (fmpz_is_strong_probabprime(p, b))
           count++;

        fmpz_clear(p);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(F);
    }

    result = (count < flint_test_multiplier());
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("count = %ld\n", count);
        fflush(stdout);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
