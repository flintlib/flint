/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_set_trunc, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        slong n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);

        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 50);

        fmpq_poly_set_trunc(b, a, n);
        fmpq_poly_set(c, a);
        fmpq_poly_truncate(c, n);

        result = (fmpq_poly_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_print(a), flint_printf("\n\n");
            fmpq_poly_print(b), flint_printf("\n\n");
            fmpq_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_set_trunc(a, a, n);

        result = (fmpq_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n");
            fmpq_poly_print(a), flint_printf("\n\n");
            fmpq_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
