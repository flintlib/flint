/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_divexact, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, q;
        int aliasing;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, 1 + n_randint(state, 100), 1 + n_randint(state, 200));
        fmpz_poly_randtest_not_zero(b, state, 1 + n_randint(state, 100), 1 + n_randint(state, 200));
        fmpz_poly_mul(c, a, b);
        aliasing = n_randint(state, 3);

        if (aliasing == 0)
        {
            fmpz_poly_divexact(q, c, b);
        }
        else if (aliasing == 1)
        {
            fmpz_poly_set(q, c);
            fmpz_poly_divexact(q, q, b);
        }
        else
        {
            fmpz_poly_set(q, b);
            fmpz_poly_divexact(q, c, q);
        }

        result = fmpz_poly_equal(q, a);

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
