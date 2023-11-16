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

TEST_FUNCTION_START(qqbar_div, state)
{
    slong iter;

    /* Check division with degree-1 terms, large coefficients */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        do {
            qqbar_randtest(x, state, 20, 100);
        } while (qqbar_is_zero(x));
        do {
            qqbar_randtest(y, state, 1, 100);
        } while (qqbar_is_zero(y));
        qqbar_randtest(z, state, 1, 100);

        /* check z / (x / y) = (z / x) * y */
        qqbar_div(a, x, y);
        qqbar_div(a, z, a);
        qqbar_div(b, z, x);
        qqbar_mul(b, b, y);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    /* Check division with degree-1 terms, small coefficients */
    for (iter = 0; iter < 100; iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        do {
            qqbar_randtest(x, state, 30, 10);
        } while (qqbar_is_zero(x));
        do {
            qqbar_randtest(y, state, 1, 10);
        } while (qqbar_is_zero(y));
        qqbar_randtest(z, state, 1, 10);

        /* check z / (x / y) = (z / x) * y */
        qqbar_div(a, x, y);
        qqbar_div(a, z, a);
        qqbar_div(b, z, x);
        qqbar_mul(b, b, y);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    /* Check division with higher-degree terms */
    for (iter = 0; iter < 100; iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        do {
            qqbar_randtest(x, state, 6, 10);
        } while (qqbar_is_zero(x));
        do {
            qqbar_randtest(y, state, 6, 10);
        } while (qqbar_is_zero(y));
        qqbar_randtest(z, state, 2, 10);

        /* check z / (x / y) = (z / x) * y */
        qqbar_div(a, x, y);
        qqbar_div(a, z, a);
        qqbar_div(b, z, x);
        qqbar_mul(b, b, y);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    /* Check roots of rational numbers, which are special-cased. */
    for (iter = 0; iter < 100; iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        do {
            qqbar_randtest(x, state, 1, 10);
        } while (qqbar_is_zero(x));
        qqbar_pow_ui(x, x, 1 + n_randint(state, 3));
        qqbar_abs(x, x);
        qqbar_root_ui(x, x, 1 + n_randint(state, 10));

        do {
            qqbar_randtest(y, state, 1, 10);
        } while (qqbar_is_zero(y));
        qqbar_pow_ui(y, y, 1 + n_randint(state, 3));
        qqbar_abs(y, y);
        qqbar_root_ui(y, y, 1 + n_randint(state, 10));

        qqbar_randtest(z, state, 2, 10);

        /* check z / (x / y) = (z / x) * y */
        qqbar_div(a, x, y);
        qqbar_div(a, z, a);
        qqbar_div(b, z, x);
        qqbar_mul(b, b, y);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    TEST_FUNCTION_END(state);
}
