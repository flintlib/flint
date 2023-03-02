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

    flint_printf("log_ui_from_prev....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t z1, z2, z3;
        ulong n1, n2;
        slong prec, accuracy;

        prec = 2 + n_randint(state, 3000);
        n1 = n_randint(state, 100000);
        n2 = n1 + 1 + n_randint(state, 1000);

        arb_init(z1);
        arb_init(z2);
        arb_init(z3);

        arb_log_ui(z1, n1, prec);
        arb_log_ui(z2, n2, prec);

        arb_log_ui_from_prev(z3, n2, z1, n1, prec);

        if (!arb_overlaps(z2, z3))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd, n1 = %wu, n2 = %wu\n\n", prec, n1, n2);
            flint_printf("z1 = "); arb_printd(z1, prec / 3.33); flint_printf("\n\n");
            flint_printf("z2 = "); arb_printd(z2, prec / 3.33); flint_printf("\n\n");
            flint_printf("z3 = "); arb_printd(z3, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        accuracy = arb_rel_accuracy_bits(z3);

        if (accuracy < prec - 4)
        {
            flint_printf("FAIL: accuracy = %wd, prec = %wd\n\n", accuracy, prec);
            flint_printf("n1 = %wu, n2 = %wu\n\n", n1, n2);
            flint_printf("z3 = "); arb_printd(z3, prec / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(z1);
        arb_clear(z2);
        arb_clear(z3);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
