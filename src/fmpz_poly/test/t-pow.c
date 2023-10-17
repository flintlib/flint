/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_pow, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        ulong exp;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(b, state, n_randint(state, 10), 100);

        exp = n_randtest(state) % UWORD(20);

        fmpz_poly_pow(a, b, exp);
        fmpz_poly_pow(b, b, exp);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp = %wu\n", exp);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Compare with repeated multiplications by the case */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        ulong exp;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 10), 100);

        exp = n_randtest(state) % UWORD(20);

        fmpz_poly_pow(a, b, exp);

        if (exp == UWORD(0))
        {
            fmpz_poly_set_ui(c, 1);
        }
        else
        {
            ulong j;
            fmpz_poly_set(c, b);

            for (j = 1; j < exp; j++)
                fmpz_poly_mul(c, c, b);
        }

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp = %wu\n", exp);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
