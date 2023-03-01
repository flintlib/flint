/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("erf....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec1, prec2, prec3, prec4;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);
        prec3 = 2 + n_randint(state, 1000);
        prec4 = 2 + n_randint(state, 1000);

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest_special(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_special(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_special(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_erf_asymp(b, a, 0, prec1, prec3);
                break;
            case 1:
                acb_hypgeom_erf_1f1a(b, a, prec1);
                break;
            case 2:
                acb_hypgeom_erf_1f1b(b, a, prec1);
                break;
            default:
                acb_hypgeom_erf(b, a, prec1);
        }

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_erf_asymp(c, a, 0, prec2, prec4);
                break;
            case 1:
                acb_hypgeom_erf_1f1a(c, a, prec2);
                break;
            case 2:
                acb_hypgeom_erf_1f1b(c, a, prec2);
                break;
            default:
                acb_hypgeom_erf(c, a, prec2);
        }

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
