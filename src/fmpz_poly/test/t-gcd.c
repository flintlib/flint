/*
    Copyright (C) 2009 William Hart
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
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_gcd, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd(a, b, c);
        fmpz_poly_gcd(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and b):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_gcd(a, b, c);
        fmpz_poly_gcd(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and c):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that a divides GCD(af, ag) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, d, f, g, q, r;

        fmpz_poly_init(a);
        fmpz_poly_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 24) + 1, 24);
        fmpz_poly_randtest(f, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);

        fmpz_poly_mul(f, a, f);
        fmpz_poly_mul(g, a, g);
        fmpz_poly_gcd(d, f, g);

        fmpz_poly_divrem_divconquer(q, r, d, a);

        result = (r->length == WORD(0));
        if (!result)
        {
            flint_printf("FAIL (check a | gcd(af, ag)):\n");
            fmpz_poly_print(f), flint_printf("\n");
            fmpz_poly_print(g), flint_printf("\n");
            fmpz_poly_print(a), flint_printf("\n");
            fmpz_poly_print(d), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}
