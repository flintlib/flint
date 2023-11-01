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

TEST_FUNCTION_START(fmpz_poly_equal_trunc, state)
{
    int i, result;

    /* equal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        fmpz_t c;
        slong n, j;

        fmpz_init(c);
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 100);

        for (j = 0; j < n; j++)
        {
           fmpz_poly_get_coeff_fmpz(c, a, j);
           fmpz_poly_set_coeff_fmpz(b, j, c);
        }

        result = (fmpz_poly_equal_trunc(a, b, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* unequal polynomials */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        fmpz_t c;
        slong m, n, j;

        fmpz_init(c);
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 100) + 1;
        m = n_randint(state, n);

        for (j = 0; j < n; j++)
        {
           fmpz_poly_get_coeff_fmpz(c, a, j);
           fmpz_poly_set_coeff_fmpz(b, j, c);
        }
        fmpz_poly_get_coeff_fmpz(c, b, m);
        fmpz_add_ui(c, c, 1);
        fmpz_poly_set_coeff_fmpz(b, m, c);

        result = (!fmpz_poly_equal_trunc(a, b, n));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            flint_printf("m = %wd\n", m);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
