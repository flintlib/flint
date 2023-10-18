/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_inflate, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t poly1, poly2, poly3, xp;
        fmpz_t one;
        ulong inflation;

        fmpz_poly_init(poly1);
        fmpz_poly_init(poly2);
        fmpz_poly_init(poly3);
        fmpz_poly_init(xp);

        fmpz_poly_randtest(poly1, state, n_randint(state, 20), n_randint(state, 200));
        inflation = n_randint(state, 10);

        fmpz_poly_inflate(poly2, poly1, inflation);

        fmpz_init(one);
        fmpz_one(one);
        fmpz_poly_set_coeff_fmpz(xp, inflation, one);
        fmpz_poly_compose(poly3, poly1, xp);
        fmpz_clear(one);

        if (!fmpz_poly_equal(poly2, poly3))
        {
            flint_printf("FAIL: not equal to compose (inflation = %wu)\n", inflation);
            flint_printf("poly1:\n");
            fmpz_poly_print(poly1);
            flint_printf("\n\n");
            flint_printf("poly2:\n");
            fmpz_poly_print(poly2);
            flint_printf("\n\n");
            flint_printf("poly3:\n");
            fmpz_poly_print(poly3);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_inflate(poly1, poly1, inflation);
        if (!fmpz_poly_equal(poly1, poly2))
        {
            flint_printf("FAIL: aliasing (inflation = %wu)\n", inflation);
            flint_printf("poly1:\n");
            fmpz_poly_print(poly1);
            flint_printf("\n\n");
            flint_printf("poly2:\n");
            fmpz_poly_print(poly2);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(poly1);
        fmpz_poly_clear(poly2);
        fmpz_poly_clear(poly3);
        fmpz_poly_clear(xp);
    }

    TEST_FUNCTION_END(state);
}
