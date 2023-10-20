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
#include "long_extras.h"
#include "fmpz_poly_q.h"

TEST_FUNCTION_START(fmpz_poly_q_scalar_mul_si, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a, b;
        slong x;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        x = z_randtest(state);

        fmpz_poly_q_scalar_mul_si(a, b, x);
        fmpz_poly_q_scalar_mul_si(b, b, x);

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
    }

    /* Check that x (a + b) == x * a + x * b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_q_t a, b, c, d;
        slong x;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        x = z_randtest(state);

        fmpz_poly_q_scalar_mul_si(c, a, x);
        fmpz_poly_q_scalar_mul_si(d, b, x);
        fmpz_poly_q_add(d, c, d);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_scalar_mul_si(c, c, x);

        result = fmpz_poly_q_equal(c, d) && fmpz_poly_q_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_q_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_q_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_q_print(c), flint_printf("\n\n");
            flint_printf("d = "), fmpz_poly_q_print(d), flint_printf("\n\n");
            flint_printf("x = %wd\n\n", x);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
    }

    TEST_FUNCTION_END(state);
}
