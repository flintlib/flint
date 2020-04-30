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

    flint_printf("pow_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* Check x^m x^n = x^(m+n) */
    for (iter = 0; iter < 100; iter++)
    {
        ca_qqbar_t x, xm, xn, xmxn, xmn;
        ulong m, n;

        ca_qqbar_init(x);
        ca_qqbar_init(xm);
        ca_qqbar_init(xn);
        ca_qqbar_init(xmxn);
        ca_qqbar_init(xmn);

        ca_qqbar_randtest(x, state, 4, 10);
        m = n_randint(state, 7);
        n = n_randint(state, 7);

        ca_qqbar_pow_ui(xm, x, m);
        ca_qqbar_pow_ui(xn, x, n);
        ca_qqbar_mul(xmxn, xm, xn);

        ca_qqbar_pow_ui(xmn, x, m + n);

        if (!ca_qqbar_equal(xmxn, xmn))
        {
            flint_printf("FAIL!\n");
            flint_printf("m = %wu\n\n", m);
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("xm = "); ca_qqbar_print(xm); flint_printf("\n\n");
            flint_printf("xn = "); ca_qqbar_print(xn); flint_printf("\n\n");
            flint_printf("xmxn = "); ca_qqbar_print(xmxn); flint_printf("\n\n");
            flint_printf("xmn = "); ca_qqbar_print(xmn); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(xm);
        ca_qqbar_clear(xn);
        ca_qqbar_clear(xmxn);
        ca_qqbar_clear(xmn);
    }

    /* Check (xy)^n = x^n y^n */
    for (iter = 0; iter < 100; iter++)
    {
        ca_qqbar_t x, y, xyn, xn, yn, xnyn;
        ulong n;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(xyn);
        ca_qqbar_init(xn);
        ca_qqbar_init(yn);
        ca_qqbar_init(xnyn);

        ca_qqbar_randtest(x, state, 4, 10);
        ca_qqbar_randtest(y, state, 4, 10);
        n = n_randint(state, 7);

        ca_qqbar_mul(xyn, x, y);
        ca_qqbar_pow_ui(xyn, xyn, n);

        ca_qqbar_pow_ui(xn, x, n);
        ca_qqbar_pow_ui(yn, y, n);
        ca_qqbar_mul(xnyn, xn, yn);

        if (!ca_qqbar_equal(xyn, xnyn))
        {
            flint_printf("FAIL!\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("xyn = "); ca_qqbar_print(xyn); flint_printf("\n\n");
            flint_printf("xn = "); ca_qqbar_print(xn); flint_printf("\n\n");
            flint_printf("yn = "); ca_qqbar_print(yn); flint_printf("\n\n");
            flint_printf("xnyn = "); ca_qqbar_print(xnyn); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(xyn);
        ca_qqbar_clear(xn);
        ca_qqbar_clear(yn);
        ca_qqbar_clear(xnyn);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

