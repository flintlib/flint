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

TEST_FUNCTION_START(fmpq_poly_tan_series, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_canonicalise(a);

        fmpq_poly_tan_series(b, a, n);
        fmpq_poly_tan_series(a, a, n);

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

    /* Check atan(tan(a)) = a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, tana, atantana;  /* but what is an atantana? */
        slong n = n_randint(state, 80) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(tana);
        fmpq_poly_init(atantana);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 60) + 1, 80);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_tan_series(tana, a, n);
        fmpq_poly_atan_series(atantana, tana, n);
        fmpq_poly_truncate(a, n);

        cflags |= fmpq_poly_is_canonical(tana) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(atantana) ? 0 : 2;
        result = (fmpq_poly_equal(atantana, a) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("tan(a) = "), fmpq_poly_debug(tana), flint_printf("\n\n");
            flint_printf("atan(tan(a)) = "), fmpq_poly_debug(atantana), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(tana);
        fmpq_poly_clear(atantana);
    }

    TEST_FUNCTION_END(state);
}
