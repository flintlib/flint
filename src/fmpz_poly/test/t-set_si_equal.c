/*
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
#include "long_extras.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_set_si_equal, state)
{
    int i, result;

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        slong n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        n = z_randtest(state);

        fmpz_poly_set_si(a, n);
        fmpz_poly_set(b, a);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n\n", n);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        slong m, n;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        m = z_randtest(state);
        n = z_randtest(state);
        while (m == n)
            n = z_randtest(state);
        fmpz_poly_set_si(a, m);
        fmpz_poly_set_si(b, n);

        result = (!fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("m = %wd\n\n", m);
            flint_printf("n = %wd\n\n", n);
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
