/*
    Copyright (C) 2021 Mathieu Gouttenoire

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_nth_derivative, state)
{
    ulong nth, cflags = UWORD(0);
    int i, j, result;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        nth = n_randint(state, 100);

        fmpq_poly_nth_derivative(b, a, nth);
        fmpq_poly_nth_derivative(a, a, nth);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check if derivative is correct */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        nth = n_randint(state, 100);

        fmpq_poly_nth_derivative(b, a, nth);
        for (j = 0; j < nth; j ++)
        {
            fmpq_poly_derivative(a, a);
        }

        result = (fmpq_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_print(a), flint_printf("\n\n");
            fmpq_poly_print(b), flint_printf("\n\n");
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
