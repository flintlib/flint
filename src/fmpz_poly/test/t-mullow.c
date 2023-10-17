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

TEST_FUNCTION_START(fmpz_poly_mullow, state)
{
    int i, result;

    /* Compare with truncated product of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        slong trunc;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        trunc = n_randint(state, 50);
        fmpz_poly_randtest(b, state, trunc, 200);
        fmpz_poly_randtest(c, state, trunc, 200);

        fmpz_poly_mullow(a, b, c, trunc);
        fmpz_poly_mul(b, b, c);
        fmpz_poly_truncate(b, trunc);

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
        fmpz_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;
        slong trunc;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 500);
        fmpz_poly_set(c, b);

        trunc = n_randint(state, 50);
        fmpz_poly_mullow(a, b, c, trunc);
        fmpz_poly_mullow(d, b, b, trunc);

        result = (fmpz_poly_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
