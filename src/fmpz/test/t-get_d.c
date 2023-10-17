/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_get_d, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, y;
        mpz_t z;
        double a, b;

        fmpz_init(x);
        fmpz_init(y);
        mpz_init(z);

        fmpz_randtest(x, state, 200);
        fmpz_get_mpz(z, x);

        a = fmpz_get_d(x);
        b = mpz_get_d(z);

        result = (a == b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("a = %f\n", a);
            flint_printf("b = %f\n", b);
            fflush(stdout);
            flint_abort();
        }

        a = a * (n_randtest(state) / (double) n_randtest_not_zero(state));

        fmpz_set_d(x, a);
        mpz_set_d(z, a);

        fmpz_set_mpz(y, z);
        result = fmpz_equal(x, y);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            flint_printf("y = "), fmpz_print(y), flint_printf("\n");
            flint_printf("a = %f\n", a);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(y);
        mpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}
