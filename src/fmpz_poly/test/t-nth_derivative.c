/*
    Copyright (C) 2021 Mathieu Gouttenoire

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

TEST_FUNCTION_START(fmpz_poly_nth_derivative, state)
{
    ulong nth;
    int i, j, result;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);

        nth = n_randint(state, 100);

        fmpz_poly_nth_derivative(b, a, nth);
        fmpz_poly_nth_derivative(a, a, nth);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Check if derivative is correct */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);

        nth = n_randint(state, 100);

        fmpz_poly_nth_derivative(b, a, nth);
        for (j = 0; j < nth; j ++)
        {
            fmpz_poly_derivative(a, a);
        }

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
