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

TEST_FUNCTION_START(fmpz_divisible, state)
{
    int i, result;

    /* Compare with GMP:  random */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;
        int e, f;

        fmpz_init(a);
        fmpz_init(b);
        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(c, a);
        fmpz_get_mpz(d, b);

        e = fmpz_divisible(b, a);
        f = mpz_divisible_p(d, c);

        result = (e == f);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
    }

    /* Compare with GMP:  b a multiple of a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t c, d;
        int e, f;

        fmpz_init(a);
        fmpz_init(b);
        mpz_init(c);
        mpz_init(d);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        fmpz_mul(b, a, b);

        fmpz_get_mpz(c, a);
        fmpz_get_mpz(d, b);

        e = fmpz_divisible(b, a);
        f = mpz_divisible_p(d, c);

        result = (e == f && e == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("c = %Zd, d = %Zd\n", c, d);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        int b;

        fmpz_init(a);

        fmpz_randtest(a, state, 200);

        b = fmpz_divisible(a, a);

        result = (b == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_print(a);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}
