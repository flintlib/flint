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

    flint_printf("sgn_re....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, z;
        int s1, s2, s3, s4;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(z);

        ca_qqbar_randtest(x, state, 4, 10);
        ca_qqbar_i(y);
        ca_qqbar_mul(y, y, x);

        s1 = ca_qqbar_sgn_re(x);
        s2 = ca_qqbar_sgn_im(y);

        if (s1 != s2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("s1 = %d", s1); flint_printf("\n\n");
            flint_printf("s2 = %d", s2); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_randtest(x, state, 4, 10);
        ca_qqbar_randtest(y, state, 4, 10);
        if (n_randint(state, 2))
            ca_qqbar_add(x, x, y);
        else
            ca_qqbar_mul(x, x, y);

        ca_qqbar_i(y);
        ca_qqbar_mul(y, y, x);

        s1 = ca_qqbar_sgn_re(x);
        s2 = ca_qqbar_sgn_im(y);

        if (s1 != s2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("s1 = %d", s1); flint_printf("\n\n");
            flint_printf("s2 = %d", s2); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_randtest_real(x, state, 2, 100);
        ca_qqbar_randtest_real(y, state, 2, 100);

        s1 = ca_qqbar_sgn_re(x);
        s2 = ca_qqbar_sgn_re(y);

        ca_qqbar_i(z);
        ca_qqbar_mul(y, y, z);
        ca_qqbar_add(x, x, y);

        s3 = ca_qqbar_sgn_re(x);
        s4 = ca_qqbar_sgn_im(x);

        if (s1 != s3 || s2 != s4)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("s1 = %d", s1); flint_printf("\n\n");
            flint_printf("s2 = %d", s2); flint_printf("\n\n");
            flint_printf("s3 = %d", s3); flint_printf("\n\n");
            flint_printf("s4 = %d", s4); flint_printf("\n\n");
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

