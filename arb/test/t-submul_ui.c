/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "flint/ulong_extras.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("submul_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        ulong x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        x = n_randtest(state);

        prec = 2 + n_randint(state, 2000);

        arb_set_ui(b, x);
        arb_set(d, c);
        arb_submul_ui(c, a, x, prec);
        arb_submul(d, a, b, prec);

        if (!arb_equal(c, d))
        {
            flint_printf("FAIL\n\n");
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

    /* aliasing */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        ulong x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        x = n_randtest(state);

        prec = 2 + n_randint(state, 2000);

        arb_set_ui(b, x);
        arb_set(c, a);
        arb_submul_ui(c, a, x, prec);
        arb_submul_ui(a, a, x, prec);

        if (!arb_equal(a, c))
        {
            flint_printf("FAIL (aliasing)\n\n");
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
