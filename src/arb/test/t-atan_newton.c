/*
    Copyright (C) 2012-2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_atan_newton, state)
{
    slong iter;

    /* Check large arguments. */
    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 100);

        arb_atan_newton(b, a, prec1);
        arb_atan_newton(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* check tan(atan(x)) = x */
        arb_sin_cos(c, d, b, prec1);
        arb_div(c, c, d, prec1);

        if (!arb_contains(c, a))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    /* Higher precision + large arguments. */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 5000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_precise(a, state, 1 + n_randint(state, 5000), 100);

        arb_atan_newton(b, a, prec1);
        arb_atan_newton(c, a, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* check tan(atan(x)) = x */
        arb_sin_cos(c, d, b, prec1);
        arb_div(c, c, d, prec1);

        if (!arb_contains(c, a))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    /* Check wide arguments. */
    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_precise(a, state, 1 + n_randint(state, 1000), 100);
        arb_randtest_precise(b, state, 1 + n_randint(state, 1000), 100);
        if (n_randint(state, 2))
            arb_add(a, a, b, 2 + n_randint(state, 1000));
        arb_union(d, a, b, 2 + n_randint(state, 1000));

        arb_atan_newton(a, a, 2 + n_randint(state, 2000));
        arb_atan_newton(b, b, 2 + n_randint(state, 2000));
        arb_atan_newton(c, d, 2 + n_randint(state, 2000));

        if (!arb_overlaps(c, a) || !arb_overlaps(c, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
