/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_deflate, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t poly1, poly2, poly3;
        ulong infl1, infl, deflation;

        fmpz_poly_init(poly1);
        fmpz_poly_init(poly2);
        fmpz_poly_init(poly3);

        fmpz_poly_randtest(poly1, state, n_randint(state, 20), n_randint(state, 200));

        if (fmpz_poly_length(poly1) <= 1)
        {
            if (fmpz_poly_deflation(poly1) != fmpz_poly_length(poly1))
            {
                flint_printf("FAIL: wrong deflation for constant polynomial\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_deflate(poly2, poly1, n_randint(state, 5) + 1);
            if (!fmpz_poly_equal(poly2, poly1))
            {
                flint_printf("FAIL: constant polynomial changed on deflation\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            infl = n_randint(state, 15) + 1;
            infl1 = fmpz_poly_deflation(poly1);

            fmpz_poly_inflate(poly2, poly1, infl);

            deflation = fmpz_poly_deflation(poly2);

            if (deflation != infl * infl1)
            {
                flint_printf("FAIL: deflation = %wu, inflation: %wu, %wu\n",
                                                       deflation, infl, infl1);
                flint_printf("poly1:\n");
                fmpz_poly_print(poly1);
                flint_printf("\n\n");
                flint_printf("poly2:\n");
                fmpz_poly_print(poly2);
                flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_deflate(poly3, poly2, infl);
            if (!fmpz_poly_equal(poly3, poly1))
            {
                flint_printf("FAIL: deflation = %wu, inflation: %wu, %wu\n",
                                                       deflation, infl, infl1);
                flint_printf("Deflated polynomial not equal to input:\n");
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

            fmpz_poly_deflate(poly2, poly2, infl);
            if (!fmpz_poly_equal(poly3, poly2))
            {
                flint_printf("FAIL: aliasing\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_poly_clear(poly1);
        fmpz_poly_clear(poly2);
        fmpz_poly_clear(poly3);
    }

    TEST_FUNCTION_END(state);
}
