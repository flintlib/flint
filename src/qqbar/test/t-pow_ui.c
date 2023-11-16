/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_pow_ui, state)
{
    slong iter;

    /* Check x^m x^n = x^(m+n) */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, xm, xn, xmxn, xmn;
        ulong m, n;

        qqbar_init(x);
        qqbar_init(xm);
        qqbar_init(xn);
        qqbar_init(xmxn);
        qqbar_init(xmn);

        if (n_randint(state, 2))
        {
            qqbar_randtest(x, state, 4, 10);
        }
        else
        {
            /* check powers of rationals which are special-cased */
            qqbar_randtest(x, state, 1, 10);
            qqbar_pow_ui(x, x, 1 + n_randint(state, 3));
            qqbar_abs(x, x);
            qqbar_root_ui(x, x, 1 + n_randint(state, 10));
        }

        m = n_randint(state, 7);
        n = n_randint(state, 7);

        qqbar_pow_ui(xm, x, m);
        qqbar_pow_ui(xn, x, n);
        qqbar_mul(xmxn, xm, xn);

        qqbar_pow_ui(xmn, x, m + n);

        if (!qqbar_equal(xmxn, xmn))
        {
            flint_printf("FAIL!\n");
            flint_printf("m = %wu\n\n", m);
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("xm = "); qqbar_print(xm); flint_printf("\n\n");
            flint_printf("xn = "); qqbar_print(xn); flint_printf("\n\n");
            flint_printf("xmxn = "); qqbar_print(xmxn); flint_printf("\n\n");
            flint_printf("xmn = "); qqbar_print(xmn); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(xm);
        qqbar_clear(xn);
        qqbar_clear(xmxn);
        qqbar_clear(xmn);
    }

    /* Check (xy)^n = x^n y^n */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, xyn, xn, yn, xnyn;
        ulong n;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(xyn);
        qqbar_init(xn);
        qqbar_init(yn);
        qqbar_init(xnyn);

        if (n_randint(state, 2))
        {
            qqbar_randtest(x, state, 4, 10);
        }
        else
        {
            /* check powers of rationals which are special-cased */
            qqbar_randtest(x, state, 1, 10);
            qqbar_pow_ui(x, x, 1 + n_randint(state, 3));
            qqbar_abs(x, x);
            qqbar_root_ui(x, x, 1 + n_randint(state, 10));
        }

        qqbar_randtest(y, state, 4, 10);
        n = n_randint(state, 7);

        qqbar_mul(xyn, x, y);
        qqbar_pow_ui(xyn, xyn, n);

        qqbar_pow_ui(xn, x, n);
        qqbar_pow_ui(yn, y, n);
        qqbar_mul(xnyn, xn, yn);

        if (!qqbar_equal(xyn, xnyn))
        {
            flint_printf("FAIL!\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("xyn = "); qqbar_print(xyn); flint_printf("\n\n");
            flint_printf("xn = "); qqbar_print(xn); flint_printf("\n\n");
            flint_printf("yn = "); qqbar_print(yn); flint_printf("\n\n");
            flint_printf("xnyn = "); qqbar_print(xnyn); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(xyn);
        qqbar_clear(xn);
        qqbar_clear(yn);
        qqbar_clear(xnyn);
    }

    TEST_FUNCTION_END(state);
}
