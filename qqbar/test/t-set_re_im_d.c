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

    flint_printf("set_re_im_d...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * calcium_test_multiplier(); iter++)
    {
        double x, y;
        qqbar_t z;
        acb_t a, b;
        int ok;

        qqbar_init(z);
        acb_init(a);
        acb_init(b);

        x = d_randtest_special(state, -1100, 1100);
        y = d_randtest_special(state, -1100, 1100);

        ok =  qqbar_set_re_im_d(z, x, y);

        qqbar_get_acb(a, z, 53);

        acb_set_d_d(b, x, y);

        if (ok)
        {
            if (!acb_equal(a, b))
            {
                flint_printf("FAIL!\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_printf("a = "); acb_print(a); flint_printf("\n\n");
                flint_printf("b = "); acb_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            if (acb_is_finite(b))
            {
                flint_printf("FAIL!\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_printf("a = "); acb_print(a); flint_printf("\n\n");
                flint_printf("b = "); acb_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        qqbar_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

