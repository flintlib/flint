/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_pow, state)
{
    slong iter;

    /* check large arguments */
    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, e, f;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(e);
        acb_init(f);

        acb_randtest(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));
        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(a));

        acb_randtest(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 200));
        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(b));

        acb_pow(c, a, b, prec1);
        acb_pow(d, a, b, prec2);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_randtest(c, state, 1 + n_randint(state, 1000), 200);

        /* check a^(b+c) = a^b*a^c */
        acb_add(e, b, c, prec1);
        acb_pow(d, a, e, prec1);

        acb_pow(e, a, b, prec1);
        acb_pow(f, a, c, prec1);
        acb_mul(e, e, f, prec1);

        if (!acb_overlaps(d, e))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_printf("e = "); acb_print(e); flint_printf("\n\n");
            flint_abort();
        }

        acb_pow(c, a, b, prec1);
        acb_set(d, a);
        acb_pow(d, d, b, prec2);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: aliasing 1\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_set(d, b);
        acb_pow(d, a, d, prec2);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: aliasing 2\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
        acb_clear(f);
    }

    TEST_FUNCTION_END(state);
}
