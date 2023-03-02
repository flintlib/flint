/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_elliptic.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("rf....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t x, y, z, r1, r2;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        acb_init(x);
        acb_init(y);
        acb_init(z);
        acb_init(r1);
        acb_init(r2);

        acb_randtest_special(x, state, 1 + n_randint(state, 300), 1 + n_randint(state, 100));
        acb_randtest_special(y, state, 1 + n_randint(state, 300), 1 + n_randint(state, 100));
        acb_randtest_special(z, state, 1 + n_randint(state, 300), 1 + n_randint(state, 100));

        acb_elliptic_rf(r1, x, y, z, 0, prec1);

        switch (n_randint(state, 6))
        {
            case 0:
                acb_elliptic_rf(r2, x, y, z, 0, prec2);
                break;
            case 1:
                acb_elliptic_rf(r2, x, z, y, 0, prec2);
                break;
            case 2:
                acb_elliptic_rf(r2, y, x, z, 0, prec2);
                break;
            case 3:
                acb_elliptic_rf(r2, y, z, x, 0, prec2);
                break;
            case 4:
                acb_elliptic_rf(r2, z, x, y, 0, prec2);
                break;
            default:
                acb_elliptic_rf(r2, z, y, x, 0, prec2);
                break;
        }

        if (!acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); acb_printd(y, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        acb_clear(r1);
        acb_clear(r2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

