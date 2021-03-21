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

    flint_printf("get_fexpr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        fexpr_t e;

        qqbar_init(x);
        qqbar_init(y);
        fexpr_init(e);

        qqbar_randtest(x, state, 5, 10);

        if (qqbar_get_fexpr_formula(e, x, 0))
        {
            qqbar_set_fexpr(y, e);

            if (!qqbar_equal(x, y))
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
                flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fexpr_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

