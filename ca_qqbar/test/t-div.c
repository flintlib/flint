/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("div....");
    fflush(stdout);

    flint_randinit(state);

    /* Check division with degree-1 terms, large coefficients */
    for (iter = 0; iter < 100; iter++)
    {
        ca_qqbar_t x, y, z, a, b;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(z);
        ca_qqbar_init(a);
        ca_qqbar_init(b);

        do {
            ca_qqbar_randtest(x, state, 20, 100);
        } while (ca_qqbar_is_zero(x));
        do {
            ca_qqbar_randtest(y, state, 1, 100);
        } while (ca_qqbar_is_zero(y));
        ca_qqbar_randtest(z, state, 1, 100);

        /* check z / (x / y) = (z / x) * y */
        ca_qqbar_div(a, x, y);
        ca_qqbar_div(a, z, a);
        ca_qqbar_div(b, z, x);
        ca_qqbar_mul(b, b, y);

        if (!ca_qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); ca_qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); ca_qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); ca_qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(z);
        ca_qqbar_clear(a);
        ca_qqbar_clear(b);
    }

    /* Check division with degree-1 terms, small coefficients */
    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, z, a, b;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(z);
        ca_qqbar_init(a);
        ca_qqbar_init(b);

        do {
            ca_qqbar_randtest(x, state, 30, 10);
        } while (ca_qqbar_is_zero(x));
        do {
            ca_qqbar_randtest(y, state, 1, 10);
        } while (ca_qqbar_is_zero(y));
        ca_qqbar_randtest(z, state, 1, 10);

        /* check z / (x / y) = (z / x) * y */
        ca_qqbar_div(a, x, y);
        ca_qqbar_div(a, z, a);
        ca_qqbar_div(b, z, x);
        ca_qqbar_mul(b, b, y);

        if (!ca_qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); ca_qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); ca_qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); ca_qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(z);
        ca_qqbar_clear(a);
        ca_qqbar_clear(b);
    }

    /* Check division with higher-degree terms */
    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, z, a, b;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(z);
        ca_qqbar_init(a);
        ca_qqbar_init(b);

        do {
            ca_qqbar_randtest(x, state, 6, 10);
        } while (ca_qqbar_is_zero(x));
        do {
            ca_qqbar_randtest(y, state, 6, 10);
        } while (ca_qqbar_is_zero(y));
        ca_qqbar_randtest(z, state, 2, 10);

        /* check z / (x / y) = (z / x) * y */
        ca_qqbar_div(a, x, y);
        ca_qqbar_div(a, z, a);
        ca_qqbar_div(b, z, x);
        ca_qqbar_mul(b, b, y);

        if (!ca_qqbar_equal(a, b))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); ca_qqbar_print(z); flint_printf("\n\n");
            flint_printf("a = "); ca_qqbar_print(a); flint_printf("\n\n");
            flint_printf("b = "); ca_qqbar_print(b); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(z);
        ca_qqbar_clear(a);
        ca_qqbar_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

