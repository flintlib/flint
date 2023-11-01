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

TEST_FUNCTION_START(fmpq_poly_sinh_cosh_series, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c, d;
        int which;
        slong n = n_randint(state, 50) + 1;
        which = n_randint(state, 2);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(d);

        fmpq_poly_randtest(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, 0);

        if (which)
        {
            fmpq_poly_sinh_cosh_series(c, b, a, n);
            fmpq_poly_sinh_cosh_series(d, a, a, n);
        }
        else
        {
            fmpq_poly_sinh_cosh_series(b, c, a, n);
            fmpq_poly_sinh_cosh_series(a, d, a, n);
        }

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 4;
        cflags |= fmpq_poly_is_canonical(d) ? 0 : 8;
        result = (fmpq_poly_equal(a, b) && fmpq_poly_equal(c, d) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasinhg %d):\n", which);
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(d);
    }

    /* Compare with sinh and cosh */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, sinha1, cosha1, sinha2, cosha2;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(sinha1);
        fmpq_poly_init(sinha2);
        fmpq_poly_init(cosha1);
        fmpq_poly_init(cosha2);

        fmpq_poly_randtest(a, state, n_randint(state, 60) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, 0);

        fmpq_poly_sinh_cosh_series(sinha1, cosha1, a, n);

        fmpq_poly_sinh_series(sinha2, a, n);
        fmpq_poly_cosh_series(cosha2, a, n);

        cflags |= fmpq_poly_is_canonical(sinha1) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(cosha1) ? 0 : 2;
        result = (fmpq_poly_equal(sinha1, sinha2) && fmpq_poly_equal(cosha1, cosha2) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("sinh(a) = "), fmpq_poly_debug(sinha1), flint_printf("\n\n");
            flint_printf("cosh(a) = "), fmpq_poly_debug(cosha1), flint_printf("\n\n");
            flint_printf("sinh(a) = "), fmpq_poly_debug(sinha2), flint_printf("\n\n");
            flint_printf("cosh(a) = "), fmpq_poly_debug(cosha2), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(sinha1);
        fmpq_poly_clear(sinha2);
        fmpq_poly_clear(cosha1);
        fmpq_poly_clear(cosha2);
    }

    TEST_FUNCTION_END(state);
}
