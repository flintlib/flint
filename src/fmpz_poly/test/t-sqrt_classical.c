/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2018 William Hart

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

TEST_FUNCTION_START(fmpz_poly_sqrt_classical, state)
{
    int i;

    /* Test aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b;
        int square1, square2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);

        fmpz_poly_randtest(a, state, 1 + n_randint(state, 20),
            1 + n_randint(state, 200));

        if (n_randint(state, 2))
            fmpz_poly_sqr(a, a);

        square1 = fmpz_poly_sqrt_classical(b, a);
        square2 = fmpz_poly_sqrt_classical(a, a);

        if ((square1 != square2) || (square1 && !fmpz_poly_equal(a, b)))
        {
            flint_printf("FAIL: aliasing:\n");
            flint_printf("square1 = %d, square2 = %d\n\n", square1, square2);
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        int square;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);

        fmpz_poly_randtest(a, state, 1 + n_randint(state, 20),
            1 + n_randint(state, 200));
        fmpz_poly_sqr(b, a);
        square = fmpz_poly_sqrt_classical(c, b);

        if (!square)
        {
            flint_printf("FAIL: square reported nonsquare:\n");
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_poly_is_zero(c) &&
            fmpz_sgn(fmpz_poly_get_coeff_ptr(c, fmpz_poly_degree(c))) < 0)
        {
            flint_printf("FAIL: leading coefficient not positive:\n");
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_sqr(c, c);
        if (!fmpz_poly_equal(c, b))
        {
            flint_printf("FAIL: sqrt(b)^2 != b:\n");
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Test "almost" squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;
        fmpz_t t;
        slong j;
        int square;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_init(t);

        fmpz_poly_randtest_not_zero(a, state, 1 + n_randint(state, 20),
            1 + n_randint(state, 200));
        fmpz_poly_sqr(b, a);

        j = n_randint(state, fmpz_poly_length(b));
        fmpz_randtest_not_zero(t, state, 1 + n_randint(state, 100));
        fmpz_add(b->coeffs + j, b->coeffs + j, t);
        _fmpz_poly_normalise(b);

        square = fmpz_poly_sqrt_classical(c, b);

        if (square)
        {
            fmpz_poly_sqr(c, c);
            if (!fmpz_poly_equal(c, b))
            {
                flint_printf("FAIL: sqrt(b)^2 != b:\n");
                flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
                flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
                flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(t);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
