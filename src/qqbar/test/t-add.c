/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("add....");
    fflush(stdout);

    flint_randinit(state);

    /* Check addition with degree-1 terms, large coefficients */
    for (iter = 0; iter < 100 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        qqbar_randtest(x, state, 20, 100);
        qqbar_randtest(y, state, 1, 100);
        qqbar_randtest(z, state, 1, 100);

        /* check (x + y) + z = x + (y + z) */
        qqbar_add(a, x, y);
        qqbar_add(a, a, z);
        qqbar_add(b, y, z);
        qqbar_add(b, x, b);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    /* Check addition with degree-1 terms, small coefficients */
    for (iter = 0; iter < 100 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        qqbar_randtest(x, state, 30, 10);
        qqbar_randtest(y, state, 1, 10);
        qqbar_randtest(z, state, 1, 10);

        /* check (x + y) + z = x + (y + z) */
        qqbar_add(a, x, y);
        qqbar_add(a, a, z);
        qqbar_add(b, y, z);
        qqbar_add(b, x, b);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    /* Check addition with higher-degree terms */
    for (iter = 0; iter < 100 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        qqbar_randtest(x, state, 6, 10);
        qqbar_randtest(y, state, 6, 10);
        qqbar_randtest(z, state, 2, 10);

        /* check (x + y) + z = x + (y + z) */
        qqbar_add(a, x, y);
        qqbar_add(a, a, z);
        qqbar_add(b, y, z);
        qqbar_add(b, x, b);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    /* More iterations, low degree */
    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, a, b;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(a);
        qqbar_init(b);

        qqbar_randtest(x, state, 4, 10);
        qqbar_randtest(y, state, 3, 10);
        qqbar_randtest(z, state, 2, 10);

        /* check (x + y) + z = x + (y + z) */
        qqbar_add(a, x, y);
        qqbar_add(a, a, z);
        qqbar_add(b, y, z);
        qqbar_add(b, x, b);

        if (!qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(a);
        qqbar_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

