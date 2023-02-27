/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("cmp_2exp_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        slong bits, e;
        arf_t x, y;
        int cmp1, cmp2;

        bits = 2 + n_randint(state, 1000);
        e = n_randtest(state);

        arf_init(x);
        arf_init(y);

        if (iter % 10 == 0)
            arf_set_ui_2exp_si(x, 1, e);
        else
            arf_randtest_special(x, state, bits, 100);

        arf_set_ui_2exp_si(y, 1, e);

        cmp1 = arf_cmp(x, y);
        cmp2 = arf_cmp_2exp_si(x, e);

        if (cmp1 != cmp2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

