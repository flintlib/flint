/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("tan....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t x, y, a, b, c, d;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(y);
        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);

        acb_randtest(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest(y, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        acb_tan(a, x, prec1);
        acb_tan(b, x, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        /* check tan(x+y) = (tan(x) + tan(y)) / (1 - tan(x) tan(y)) */
        acb_add(b, x, y, prec1);
        acb_tan(b, b, prec1);

        acb_tan(c, y, prec1);
        acb_add(d, a, c, prec1);
        acb_mul(c, a, c, prec1);
        acb_sub_ui(c, c, 1, prec1);
        acb_neg(c, c);
        acb_div(d, d, c, prec1);

        if (!acb_overlaps(b, d))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("y = "); acb_print(y); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_tan(x, x, prec1);

        if (!acb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

