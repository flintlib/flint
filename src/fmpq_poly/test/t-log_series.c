/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_log_series, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(1));

        fmpq_poly_canonicalise(a);

        fmpq_poly_log_series(b, a, n);
        fmpq_poly_log_series(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check log(a*b) = log(a) + log(b) */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, ab, loga, logb, logab, loga_logb;
        slong n = n_randint(state, 80) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(ab);
        fmpq_poly_init(loga);
        fmpq_poly_init(logb);
        fmpq_poly_init(logab);
        fmpq_poly_init(loga_logb);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 80) + 1, 80);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(1));

        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 80) + 1, 80);
        fmpq_poly_set_coeff_ui(b, 0, UWORD(1));

        fmpq_poly_mullow(ab, a, b, n);

        fmpq_poly_log_series(logab, ab, n);
        fmpq_poly_log_series(loga, a, n);
        fmpq_poly_log_series(logb, b, n);
        fmpq_poly_add(loga_logb, loga, logb);

        cflags |= fmpq_poly_is_canonical(loga) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(logb) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(logab) ? 0 : 4;
        result = (fmpq_poly_equal(logab, loga_logb) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("log(a) = "), fmpq_poly_debug(loga), flint_printf("\n\n");
            flint_printf("log(b) = "), fmpq_poly_debug(logb), flint_printf("\n\n");
            flint_printf("log(ab) = "), fmpq_poly_debug(logab), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(ab);
        fmpq_poly_clear(loga);
        fmpq_poly_clear(logb);
        fmpq_poly_clear(logab);
        fmpq_poly_clear(loga_logb);
    }

    TEST_FUNCTION_END(state);
}
