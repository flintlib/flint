/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_div_series, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing q and a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);

        fmpq_poly_randtest(a, state, n_randint(state, 50) + 1, 80);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 80);
        fmpq_poly_set_coeff_ui(b, 0, 1);

        fmpq_poly_div_series(q, a, b, n);
        fmpq_poly_div_series(a, a, b, n);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(a) ? 0 : 2;
        result = (fmpq_poly_equal(q, a)) && !cflags;
        if (!result)
        {
            flint_printf("FAIL (alias q and a):\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("q = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
    }

    /* Check aliasing q and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);

        fmpq_poly_randtest(a, state, n_randint(state, 50) + 1, 80);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 80);
        fmpq_poly_set_coeff_ui(b, 0, 1);

        fmpq_poly_div_series(q, a, b, n);
        fmpq_poly_div_series(b, a, b, n);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(q, b)) && !cflags;
        if (!result)
        {
            flint_printf("FAIL (alias q and b):\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("q = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
    }

    /* Check that Q * B == A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, p, q;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(p);
        fmpq_poly_init(q);

        fmpq_poly_randtest(a, state, n_randint(state, 50) + 1, 80);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 80);
        fmpq_poly_set_coeff_ui(b, 0, 1);

        fmpq_poly_div_series(q, a, b, n);
        fmpq_poly_mullow(p, q, b, n);

        fmpq_poly_truncate(a, n);

        cflags |= fmpq_poly_is_canonical(q) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(p) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(a) ? 0 : 4;
        result = (fmpq_poly_equal(p, a)) && !cflags;
        if (!result)
        {
            flint_printf("FAIL (check Q * B = A):\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("p = "), fmpq_poly_debug(p), flint_printf("\n\n");
            flint_printf("q = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(p);
        fmpq_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}
