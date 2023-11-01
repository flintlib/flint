/*
    Copyright (C) 2010, 2020 William Hart

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

TEST_FUNCTION_START(fmpz_poly_mul_SS_precache, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
	fmpz_poly_mul_precache_t pre;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 50), 200);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

        fmpz_poly_mul_SS_precache_init(pre, 50, 200, c);

        fmpz_poly_mul_SS(a, b, c);
        fmpz_poly_mul_SS_precache(b, b, pre);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

	fmpz_poly_mul_precache_clear(pre);

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Compare with mul_KS */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, d;
        fmpz_poly_mul_precache_t pre;
	int k;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);
        fmpz_poly_randtest(c, state, n_randint(state, 50), 200);

	fmpz_poly_mul_SS_precache_init(pre, 50, 200, c);

	for (k = 0; k < 3; k++)
	{
	   fmpz_poly_randtest(b, state, n_randint(state, 50), 200);

           fmpz_poly_mul_KS(a, b, c);
           fmpz_poly_mul_SS_precache(d, b, pre);

           result = (fmpz_poly_equal(a, d));
           if (!result)
           {
               flint_printf("FAIL:\n");
               fmpz_poly_print(a), flint_printf("\n\n");
               fmpz_poly_print(d), flint_printf("\n\n");
               fflush(stdout);
               flint_abort();
           }
	}

	fmpz_poly_mul_precache_clear(pre);

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
