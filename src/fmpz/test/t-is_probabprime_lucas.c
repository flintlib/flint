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

TEST_FUNCTION_START(fmpz_is_probabprime_lucas, state)
{
    int i, result, count = 0;

    /* test primes always pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, F;

        fmpz_init(p);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 330) + 2);
           fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p) || fmpz_cmp_ui(p, 2) == 0);

        result = fmpz_is_probabprime_lucas(p);
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

    /* test composites rarely pass */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p, a, F;

        fmpz_init(p);
        fmpz_init(a);
        fmpz_init(F);

        do {
           fmpz_randbits(p, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(p, 2) < 0);
        do {
           fmpz_randbits(a, state, n_randint(state, 100) + 2);
        } while (fmpz_cmp_ui(a, 2) < 0);

        fmpz_mul(p, p, a);

        if (fmpz_is_probabprime_lucas(p))
           count++;

        fmpz_clear(p);
        fmpz_clear(a);
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
