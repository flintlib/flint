/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

TEST_FUNCTION_START(fmpz_poly_derivative, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpz_poly_derivative(b, a);
        fmpz_poly_derivative(a, a);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Check constants have derivative zero */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 2), 200);

        fmpz_poly_derivative(b, a);

        result = (b->length == UWORD(0));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Check (f g)' = f' g + f g' */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d, lhs, rhs;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_init(lhs);
        fmpz_poly_init(rhs);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest(b, state, n_randint(state, 100), 200);

        fmpz_poly_mul(lhs, a, b);
        fmpz_poly_derivative(lhs, lhs);
        fmpz_poly_derivative(c, a);
        fmpz_poly_derivative(d, b);
        fmpz_poly_mul(c, c, b);
        fmpz_poly_mul(d, a, d);
        fmpz_poly_add(rhs, c, d);

        result = fmpz_poly_equal(lhs, rhs);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
        fmpz_poly_clear(lhs);
        fmpz_poly_clear(rhs);
    }

    TEST_FUNCTION_END(state);
}
