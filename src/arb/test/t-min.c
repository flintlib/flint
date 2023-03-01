/*
    Copyright (C) 2016 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* sample (x, y) so that x /in y */
void _sample_arf_in_arb(arf_t x, arb_t y, flint_rand_t state)
{
    slong bits, prec, expbits;
    arf_t a, b;
    slong i, n;

    arf_init(a);
    arf_init(b);

    bits = 2 + n_randint(state, 1000);
    prec = 2 + n_randint(state, 1000);
    expbits = n_randint(state, 14);
    n = n_randint(state, 3);

    arf_randtest(x, state, bits, expbits);
    arf_set(a, x);
    arf_set(b, x);
    for (i = 0; i < n; i++)
    {
        arf_randtest(x, state, bits, expbits);
        arf_min(a, a, x);
        arf_max(b, b, x);
    }

    arb_set_interval_arf(y, a, b, prec);

    arf_clear(a);
    arf_clear(b);
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("min....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arf_t a, b, c;
        arb_t x, y, z;
        slong prec;

        arf_init(a);
        arf_init(b);
        arf_init(c);

        arb_init(x);
        arb_init(y);
        arb_init(z);

        _sample_arf_in_arb(a, x, state);
        _sample_arf_in_arb(b, y, state);

        prec = 2 + n_randint(state, 200);

        arf_min(c, a, b);
        arb_min(z, x, y, prec);

        if (!arb_contains_arf(x, a) ||
            !arb_contains_arf(y, b) ||
            !arb_contains_arf(z, c))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("c = "); arf_print(c); flint_printf("\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        /* aliasing */
        {
            int alias;

            if (n_randint(state, 2))
            {
                arb_min(x, x, y, prec);
                alias = arb_equal(x, z);
            }
            else
            {
                arb_min(y, x, y, prec);
                alias = arb_equal(y, z);
            }

            if (!alias)
            {
                flint_printf("FAIL: aliasing\n\n");
                flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                flint_abort();
            }
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(c);

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
