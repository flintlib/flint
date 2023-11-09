/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_pow_trunc, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        ulong exp;
        slong len;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(b, state, n_randint(state, 10), 100);

        exp = (ulong) n_randtest(state) % UWORD(20);
        len = n_randtest(state) % 20;

        fmpq_poly_pow_trunc(a, b, exp, len);
        fmpq_poly_pow_trunc(b, b, exp, len);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp = %wu, len = %wd\n", exp, len);
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Compare with repeated multiplications by the base */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        ulong exp;
        slong len;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(b, state, n_randint(state, 10), 100);
        fmpq_poly_randtest(c, state, n_randint(state, 10), 100);

        exp = (ulong) n_randtest(state) % UWORD(20);
        len = n_randtest(state) % 20;

        fmpq_poly_pow_trunc(a, b, exp, len);

        if (exp == 0)
        {
            fmpq_poly_set_ui(c, 1);
        }
        else
        {
            ulong j;
            fmpq_poly_set(c, b);
            for (j = 1; j < exp; j++)
                fmpq_poly_mullow(c, c, b, len);
        }
        fmpq_poly_truncate(c, len);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        result = (fmpq_poly_equal(a, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp = %wu, len = %wd\n", exp, len);
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("c = "), fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
