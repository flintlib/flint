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

    flint_printf("mul_2exp_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, z;
        slong e;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(z);

        ca_qqbar_randtest(x, state, 6, 10);
        e = n_randint(state, 100) - (slong) 50;

        if (n_randint(state, 2))
            ca_qqbar_mul_2exp_si(y, x, e);
        else
        {
            ca_qqbar_set(y, x);
            ca_qqbar_mul_2exp_si(y, y, e);
        }

        ca_qqbar_one(z);
        ca_qqbar_mul_2exp_si(z, z, e);
        ca_qqbar_mul(z, z, x);

        if (!ca_qqbar_equal(z, y))
        {
            flint_printf("FAIL!\n");
            flint_printf("%wd\n\n", e);
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); ca_qqbar_print(z); flint_printf("\n\n");
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

