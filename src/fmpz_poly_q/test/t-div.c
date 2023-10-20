/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly_q.h"

TEST_FUNCTION_START(fmpz_poly_q_div, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a, b, c;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest_not_zero(c, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_div(a, b, c);
        fmpz_poly_q_div(b, b, c);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a), flint_printf("\n\n");
            fmpz_poly_q_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a, b, c;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest_not_zero(c, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_div(a, b, c);
        fmpz_poly_q_div(c, b, c);

        result = (fmpz_poly_q_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a), flint_printf("\n\n");
            fmpz_poly_q_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
    }

    /* Check (c/b)+(d/b) = (c+d)/b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a1, a2, b, c, d;

        fmpz_poly_q_init(a1);
        fmpz_poly_q_init(a2);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_poly_q_randtest_not_zero(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(c, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(d, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_div(a1, c, b);
        fmpz_poly_q_div(a2, d, b);
        fmpz_poly_q_add(a1, a1, a2);

        fmpz_poly_q_add(c, c, d);
        fmpz_poly_q_div(a2, c, b);

        result = fmpz_poly_q_equal(a1, a2) && fmpz_poly_q_is_canonical(a1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a1), flint_printf("\n\n");
            fmpz_poly_q_print(a2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_q_clear(a1);
        fmpz_poly_q_clear(a2);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
    }

    TEST_FUNCTION_END(state);
}
