/*
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

TEST_FUNCTION_START(fmpz_poly_pseudo_rem_cohen, state)
{
    int i, result;

    /* Compare with q*b + r = a, no aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r1, r2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r1);
        fmpz_poly_init(r2);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);

        fmpz_poly_pseudo_divrem_cohen(q, r1, a, b);
        fmpz_poly_pseudo_rem_cohen(r2, a, b);

        result = (fmpz_poly_equal(r1, r2));
        if (!result)
        {
            flint_printf("FAIL (correctness):\n");
            flint_printf("a  = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b  = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("r1 = "), fmpz_poly_print(r1), flint_printf("\n\n");
            flint_printf("r2 = "), fmpz_poly_print(r2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r1);
        fmpz_poly_clear(r2);
    }

    /* Check r and a alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, r;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);

        fmpz_poly_pseudo_rem_cohen(r, a, b);
        fmpz_poly_pseudo_rem_cohen(a, a, b);

        result = (fmpz_poly_equal(a, r));
        if (!result)
        {
            flint_printf("FAIL (aliasing r, a):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(r);
    }

    /* Check r and b alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, r;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);

        fmpz_poly_pseudo_rem_cohen(r, a, b);
        fmpz_poly_pseudo_rem_cohen(b, a, b);

        result = (fmpz_poly_equal(b, r));
        if (!result)
        {
            flint_printf("FAIL (aliasing r, b):\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}
