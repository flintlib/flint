/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gmpcompat.h"
#include "long_extras.h"

TEST_FUNCTION_START(fmpz_divisible_si, state)
{
    int i, result;

    /* Compare with GMP:  random */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong a;
        fmpz_t b;
        mpz_t d;
        int e, f;

        fmpz_init(b);
        mpz_init(d);

        a = z_randtest_not_zero(state);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, b);

        e = fmpz_divisible_si(b, a);
        f = flint_mpz_divisible_ui_p(d, FLINT_ABS(a));

        result = (e == f);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wd, b = ", a), fmpz_print(b), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(b);
        mpz_clear(d);
    }

    /* Compare with GMP:  b a multiple of a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong a;
        fmpz_t b;
        mpz_t d;
        int e, f;

        fmpz_init(b);
        mpz_init(d);

        a = z_randtest_not_zero(state);
        fmpz_randtest(b, state, 200);
        fmpz_mul_si(b, b, a);

        fmpz_get_mpz(d, b);

        e = fmpz_divisible_si(b, a);
        f = flint_mpz_divisible_ui_p(d, FLINT_ABS(a));

        result = (e == f && e == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = %wd, b = ", a), fmpz_print(b), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(b);
        mpz_clear(d);
    }

    TEST_FUNCTION_END(state);
}
