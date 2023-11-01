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

TEST_FUNCTION_START(fmpz_is_prime_pseudosquare, state)
{
    int i, result;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        int r1, r2;

        fmpz_init(p);

        fmpz_randtest_unsigned(p, state, n_randint(state, 94) + 1);

        r1 = fmpz_is_probabprime(p);
        r2 = fmpz_is_prime_pseudosquare(p);

        result = (r1 == r2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p);
            flint_printf("r1 = %d, r2 = %d\n", r1, r2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
