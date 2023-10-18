/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_divides, state)
{
    int i, result;

    /* Check that b divides a*b and that the quotient is a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        fmpz_poly_mul(p, a, b);

        result = (fmpz_poly_divides(q, p, b) && fmpz_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
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
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        fmpz_poly_mul(p, a, b);

        result = (fmpz_poly_divides(p, p, b) && fmpz_poly_equal(p, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(p), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
    }

    /* Check aliasing of q with b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        fmpz_poly_mul(p, a, b);

        result = (fmpz_poly_divides(b, p, b) && fmpz_poly_equal(b, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(p), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
    }

    /* Check when not divisible */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q, g, s;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_init(s);
        fmpz_poly_init(g);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 200);
        do {
           fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 200);
        } while (b->length < 2);
        fmpz_poly_mul(p, a, b);
        do {
           fmpz_poly_randtest_not_zero(s, state, b->length, 200);
           fmpz_poly_gcd(g, s, b);
        } while (g->length == b->length);
        fmpz_poly_add(p, p, s);

        result = (!fmpz_poly_divides(q, p, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
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
        fmpz_poly_clear(s);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
