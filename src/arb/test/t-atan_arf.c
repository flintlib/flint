/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("atan_arf....");
    fflush(stdout);

    flint_randinit(state);

    /* self-consistency test */
    for (iter = 0; iter < 5000 * arb_test_multiplier(); iter++)
    {
        arf_t x;
        arb_t y1, y2;
        slong prec1, prec2, acc1, acc2;

        prec1 = 2 + n_randint(state, 9000);
        prec2 = 2 + n_randint(state, 9000);

        arf_init(x);
        arb_init(y1);
        arb_init(y2);

        arf_randtest_special(x, state, 1 + n_randint(state, 9000), 200);
        arb_randtest_special(y1, state, 1 + n_randint(state, 9000), 200);
        arb_randtest_special(y2, state, 1 + n_randint(state, 9000), 200);

        if (n_randint(state, 2))
            arf_add_ui(x, x, 1, 2 + n_randint(state, 9000), ARF_RND_DOWN);

        arb_atan_arf(y1, x, prec1);
        arb_atan_arf(y2, x, prec2);

        if (!arb_overlaps(y1, y2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y1 = "); arb_print(y1); flint_printf("\n\n");
            flint_printf("y2 = "); arb_print(y2); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(y1);
        acc2 = arb_rel_accuracy_bits(y2);

        if (!arf_is_nan(x))
        {
            if (acc1 < prec1 - 2 || acc2 < prec2 - 2)
            {
                flint_printf("FAIL: accuracy\n\n");
                flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
                flint_printf("acc1 = %wd, acc2 = %wd\n\n", acc1, acc2);
                flint_printf("x = "); arf_printd(x, 50); flint_printf("\n\n");
                flint_printf("y1 = "); arb_printd(y1, 50); flint_printf("\n\n");
                flint_printf("y2 = "); arb_printd(y2, 50); flint_printf("\n\n");
                flint_abort();
            }
        }

        arf_clear(x);
        arb_clear(y1);
        arb_clear(y2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

