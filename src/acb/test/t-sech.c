/*
    Copyright (C) 2017 Fredrik Johansson

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

    flint_printf("sech....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t x, a, b;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 200);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(a);
        acb_init(b);

        acb_randtest_special(x, state, 1 + n_randint(state, 300), 100);
        acb_randtest_special(a, state, 1 + n_randint(state, 300), 100);
        acb_randtest_special(b, state, 1 + n_randint(state, 300), 100);

        if (n_randint(state, 2))
        {
            acb_sech(a, x, prec1);
        }
        else
        {
            acb_set(a, x);
            acb_sech(a, a, prec1);
        }

        acb_cosh(b, x, prec2);
        acb_inv(b, b, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printn(x, 20, 0); flint_printf("\n\n");
            flint_printf("a = "); acb_printn(a, 20, 0); flint_printf("\n\n");
            flint_printf("b = "); acb_printn(b, 20, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(a);
        acb_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

