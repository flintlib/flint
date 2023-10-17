/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_is_cyclotomic, state)
{
    int i;

    /* Check detection of small cyclotomics */
    for (i = 0; i < 250; i++)
    {
        fmpz_poly_t f;
        ulong n, r;

        n = i;

        fmpz_poly_init(f);
        fmpz_poly_cyclotomic(f, n);
        r = fmpz_poly_is_cyclotomic(f);

        if (r != n)
        {
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("%wu %wu\n\n", n, r);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
    }

    /* Check detection of large cyclotomics */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f;
        ulong n, r;

        n = n_randtest(state) % 10000;

        fmpz_poly_init(f);
        fmpz_poly_cyclotomic(f, n);
        r = fmpz_poly_is_cyclotomic(f);

        if (r != n)
        {
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("%wu %wu\n\n", n, r);
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
    }

    /* Check cyclotomic + random polynomial */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g;
        ulong n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_poly_cyclotomic(f, 1 + n_randint(state, 100));
        fmpz_poly_randtest(g, state, 1 + n_randint(state, 50), 1);
        fmpz_poly_add(f, f, g);

        n = fmpz_poly_is_cyclotomic(f);

        if (n != 0)
        {
            fmpz_poly_cyclotomic(g, n);

            if (!fmpz_poly_equal(f, g))
            {
                flint_printf("FAIL\n");
                fmpz_poly_print(f); flint_printf("\n\n");
                fmpz_poly_print(g); flint_printf("\n\n");
                flint_printf("%wu\n\n", n);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
