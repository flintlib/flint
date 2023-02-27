/*
    Copyright (C) 2012 Fredrik Johansson

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

    flint_printf("rising_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* compare with fmpq */
    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;
        ulong n;
        slong i;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 1000), 10);
        arb_randtest(b, state, 1 + n_randint(state, 1000), 10);
        n = n_randint(state, 80);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 10));
        arb_rising_ui(b, a, n, 2 + n_randint(state, 1000));

        fmpq_one(y);
        for (i = 0; i < n; i++)
        {
            fmpq_set_si(z, i, 1);
            fmpq_add(z, x, z);
            fmpq_mul(y, y, z);
        }

        if (!arb_contains_fmpq(b, y))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of y and x */
    for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
        arb_t x, y;
        ulong n;
        slong prec;

        arb_init(x);
        arb_init(y);

        arb_randtest(x, state, 1 + n_randint(state, 200), 10);
        arb_randtest(y, state, 1 + n_randint(state, 200), 10);
        n = n_randint(state, 100);

        prec = 2 + n_randint(state, 1000);

        arb_rising_ui(y, x, n, prec);
        arb_rising_ui(x, x, n, prec);

        if (!arb_equal(x, y))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("n = %wu\n", n);
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
