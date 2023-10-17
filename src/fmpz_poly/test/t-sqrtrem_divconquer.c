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

TEST_FUNCTION_START(fmpz_poly_sqrtrem_divconquer, state)
{
    int i;

    /* Test aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, r1, r2;
        int sqrtrem1, sqrtrem2;
        slong len;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(r1);
        fmpz_poly_init(r2);

        len = n_randint(state, 150);

        do {
            fmpz_poly_randtest(a, state, len + 1, 1 + n_randint(state, 200));
            fmpz_poly_randtest(c, state, len, 1 + n_randint(state, 200));
        } while (fmpz_poly_length(a) == fmpz_poly_length(c));

        if (fmpz_poly_length(a) < fmpz_poly_length(c))
           fmpz_poly_swap(a, c);

        if (n_randint(state, 2))
        {
           fmpz_poly_sqr(a, a);
           fmpz_poly_add(a, a, c);
        }

        sqrtrem1 = fmpz_poly_sqrtrem_divconquer(b, r1, a);
        sqrtrem2 = fmpz_poly_sqrtrem_divconquer(a, r2, a);

        if ((sqrtrem1 != sqrtrem2) ||
            (sqrtrem1 && (!fmpz_poly_equal(a, b) || !fmpz_poly_equal(r1, r2))))
        {
            flint_printf("FAIL: aliasing:\n");
            flint_printf("sqrtrem1 = %d, sqrtrem2 = %d\n\n", sqrtrem1, sqrtrem2);
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("r1: "); fmpz_poly_print(r1); flint_printf("\n\n");
            flint_printf("r2: "); fmpz_poly_print(r2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(r1);
        fmpz_poly_clear(r2);
    }

    /* Test aliasing of a and r */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, r;
        int sqrtrem1, sqrtrem2;
        slong len;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(r);

        len = n_randint(state, 150);

        do {
            fmpz_poly_randtest(a, state, len + 1, 1 + n_randint(state, 200));
            fmpz_poly_randtest(c, state, len, 1 + n_randint(state, 200));
        } while (fmpz_poly_length(a) == fmpz_poly_length(c));

        if (fmpz_poly_length(a) < fmpz_poly_length(c))
           fmpz_poly_swap(a, c);

        if (n_randint(state, 2))
        {
           fmpz_poly_sqr(a, a);
           fmpz_poly_add(a, a, c);
        }

        sqrtrem1 = fmpz_poly_sqrtrem_divconquer(b, r, a);
        sqrtrem2 = fmpz_poly_sqrtrem_divconquer(c, a, a);

        if ((sqrtrem1 != sqrtrem2) ||
            (sqrtrem1 && (!fmpz_poly_equal(a, r) || !fmpz_poly_equal(b, c))))
        {
            flint_printf("FAIL: aliasing2:\n");
            flint_printf("sqrtrem1 = %d, sqrtrem2 = %d\n\n", sqrtrem1, sqrtrem2);
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
            flint_printf("r: "); fmpz_poly_print(r); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(r);
    }

    /* Test random squares with remainder */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, r1, r2;
        int sqrtrem;
        slong len;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(r1);
        fmpz_poly_init(r2);

        len = n_randint(state, 150);

        do {
            fmpz_poly_randtest(a, state, len + 1, 1 + n_randint(state, 200));
            fmpz_poly_randtest(r1, state, len, 1 + n_randint(state, 200));
        } while (fmpz_poly_length(a) == fmpz_poly_length(r1));

        if (fmpz_poly_length(a) < fmpz_poly_length(r1))
           fmpz_poly_swap(a, r1);

        fmpz_poly_sqr(b, a);
        fmpz_poly_add(b, b, r1);

        sqrtrem = fmpz_poly_sqrtrem_divconquer(c, r2, b);

        if (!sqrtrem)
        {
            flint_printf("FAIL: sqrtrem returns false:\n");
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
            flint_printf("r1: "); fmpz_poly_print(r1); flint_printf("\n\n");
            flint_printf("r2: "); fmpz_poly_print(r2); flint_printf("\n\n");
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
            flint_printf("r1: "); fmpz_poly_print(r1); flint_printf("\n\n");
            flint_printf("r2: "); fmpz_poly_print(r2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_sqr(c, c);
        fmpz_poly_add(c, c, r2);
        if (!fmpz_poly_equal(c, b))
        {
            flint_printf("FAIL: sqrt(b)^2 + r != b:\n");
            flint_printf("a: "); fmpz_poly_print(a); flint_printf("\n\n");
            flint_printf("b: "); fmpz_poly_print(b); flint_printf("\n\n");
            flint_printf("c: "); fmpz_poly_print(c); flint_printf("\n\n");
            flint_printf("r1: "); fmpz_poly_print(r1); flint_printf("\n\n");
            flint_printf("r2: "); fmpz_poly_print(r2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(r1);
        fmpz_poly_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
