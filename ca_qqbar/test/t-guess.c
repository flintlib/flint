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

    flint_printf("guess....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200; iter++)
    {
        ca_qqbar_t x, y;
        acb_t z;
        slong prec, max_deg, max_bits;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        acb_init(z);

        if (n_randint(state, 2))
            ca_qqbar_randtest(x, state, 4, 100);
        else
            ca_qqbar_randtest(x, state, 12, 10);

        max_deg = ca_qqbar_degree(x) + n_randint(state, 2);
        max_bits = ca_qqbar_height_bits(x) + n_randint(state, 100);

        for (prec = 64; ; prec *= 2)
        {
            ca_qqbar_get_acb(z, x, prec);

            if (ca_qqbar_guess(y, z, max_deg, max_bits, 0, prec) && ca_qqbar_equal(x, y))
                break;

            if (prec > 10000)
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_abort();
            }
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        acb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

