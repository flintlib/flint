/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("doublefac_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        ulong n;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        n = n_randtest(state);
        if (n + 1 == 0 || n + 2 == 0)
            n -= 2;

        arb_doublefac_ui(a, n, prec1);
        arb_doublefac_ui(b, n + 2, prec1);
        arb_mul_ui(c, a, n + 2, prec2);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
