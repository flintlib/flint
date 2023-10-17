/*
    Copyright (C) 2009, 2019 William Hart

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

TEST_FUNCTION_START(fmpz_poly_div_series_divconquer, state)
{
    int i, result;

    /* Check aliasing q and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q;
        slong n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));
        fmpz_poly_set_coeff_si(b, 0, n_randint(state, 2) ? 1 : -1);

        fmpz_poly_div_series_divconquer(q, a, b, n);
        fmpz_poly_div_series_divconquer(a, a, b, n);

        result = (fmpz_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL (alias q and a):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("q = "), fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
    }

    /* Check aliasing q and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q;
        slong n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));
        fmpz_poly_set_coeff_si(b, 0, n_randint(state, 2) ? 1 : -1);

        fmpz_poly_div_series_divconquer(q, a, b, n);
        fmpz_poly_div_series_divconquer(b, a, b, n);

        result = (fmpz_poly_equal(q, b));
        if (!result)
        {
            flint_printf("FAIL (alias q and b):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("q = "), fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
    }

    /* Check that Q * B == A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        slong n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));
        fmpz_poly_set_coeff_si(b, 0, n_randint(state, 2) ? 1 : -1);

        fmpz_poly_div_series_divconquer(q, a, b, n);
        fmpz_poly_mullow(p, q, b, n);

        fmpz_poly_truncate(a, n);

        result = (fmpz_poly_equal(p, a));
        if (!result)
        {
            flint_printf("FAIL (check Q * B = A):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("p = "), fmpz_poly_print(p), flint_printf("\n\n");
            flint_printf("q = "), fmpz_poly_print(q), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    /* Check that (A * B)/B == A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, p, q;
        slong n = n_randint(state, 50) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(p);
        fmpz_poly_init(q);

        fmpz_poly_randtest(a, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100))
;
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));

        while (fmpz_is_zero(b->coeffs + 0))
            fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 2 + n_randint(state, 100));

        fmpz_poly_mullow(p, a, b, n);

        fmpz_poly_div_series_divconquer(q, p, b, n);

        fmpz_poly_truncate(a, n);

        result = (fmpz_poly_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL (check (A * B)/B = A):\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("p = "), fmpz_poly_print(p), flint_printf("\n\n");
            flint_printf("q = "), fmpz_poly_print(q), flint_printf("\n\n");
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
