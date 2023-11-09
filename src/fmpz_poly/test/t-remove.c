/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_remove, state)
{
    int i, result;

    /* Check that b divides a*b and that the quotient is a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        slong e1, e2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);
        while (b->length == 0 || (b->length == 1 && fmpz_is_pm1(b->coeffs + 0)))
	   fmpz_poly_randtest_not_zero(b, state, n_randint(state, 10) + 1, 100);
        while (!fmpz_poly_is_one(p))
        {
            fmpz_poly_randtest(a, state, n_randint(state, 100), 100);
	    fmpz_poly_gcd(p, a, b);
	}
	e1 = n_randint(state, 10);
	fmpz_poly_pow(p, b, e1);
	fmpz_poly_mul(p, p, a);

        e2 = fmpz_poly_remove(q, p, b);

	result = (e2 == e1 && (e2 == 0 || fmpz_poly_equal(q, a)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("e1 = %wd, e2 = %wd\n", e1, e2);
	    fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(p), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    /* Check aliasing of q with a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        slong e1, e2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);

        while (b->length == 0 || (b->length == 1 && fmpz_is_pm1(b->coeffs + 0)))
	    fmpz_poly_randtest_not_zero(b, state, n_randint(state, 10) + 1, 100);
        while (!fmpz_poly_is_one(p))
        {
            fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
            fmpz_poly_gcd(p, a, b);
        }

        e1 = n_randint(state, 10);
        fmpz_poly_pow(p, b, e1);
        fmpz_poly_mul(p, p, a);

        e1 = fmpz_poly_remove(q, p, b);
        e2 = fmpz_poly_remove(p, p, b);

        result = (e1 == e2 && fmpz_poly_equal(p, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("e1 = %wd, e2 = %wd\n", e1, e2);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(p), flint_printf("\n\n");
            fmpz_poly_print(p), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    /* Check aliasing of q with b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        slong e1, e2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);

        while (b->length == 0 || (b->length == 1 && fmpz_is_pm1(b->coeffs + 0)))
            fmpz_poly_randtest_not_zero(b, state, n_randint(state, 10) + 1, 100);
        while (!fmpz_poly_is_one(p))
        {
            fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
            fmpz_poly_gcd(p, a, b);
        }

        e1 = n_randint(state, 10);
        fmpz_poly_pow(p, b, e1);
        fmpz_poly_mul(p, p, a);

        e1 = fmpz_poly_remove(q, p, b);
        e2 = fmpz_poly_remove(b, p, b);

        result = (e1 == e2 && fmpz_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("e1 = %wd, e2 = %wd\n", e1, e2);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(p), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
