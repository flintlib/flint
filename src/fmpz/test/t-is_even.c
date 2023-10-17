/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_is_even, state)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t g;

        fmpz_init(f);
        mpz_init(g);

        fmpz_randtest(f, state, 100);
        fmpz_get_mpz(g, f);

        result = (fmpz_is_even(f) == mpz_even_p(g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            gmp_printf("g = %Zd\n", g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        mpz_clear(g);
    }

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        mpz_t g;

        fmpz_init(f);
        mpz_init(g);

        fmpz_randtest(f, state, 100);
        fmpz_get_mpz(g, f);

        result = (fmpz_is_odd(f) == mpz_odd_p(g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_print(f), flint_printf("\n");
            gmp_printf("g = %Zd\n", g);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        mpz_clear(g);
    }

    TEST_FUNCTION_END(state);
}
