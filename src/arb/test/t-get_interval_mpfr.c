/*
    Copyright (C) 2014 Fredrik Johansson

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

    flint_printf("get_interval_mpfr....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y;
        mpfr_t aa, bb;

        arb_init(x);
        arb_init(y);
        mpfr_init2(aa, 2 + n_randint(state, 200));
        mpfr_init2(bb, 2 + n_randint(state, 200));

        arb_randtest_special(x, state, 200, 10);
        arb_get_interval_mpfr(aa, bb, x);
        arb_set_interval_mpfr(y, aa, bb, 2 + n_randint(state, 200));

        if (!arb_contains(y, x))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("aa = "); mpfr_printf("%.50Rg", aa); flint_printf("\n\n");
            flint_printf("bb = "); mpfr_printf("%.50Rg", bb); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        mpfr_clear(aa);
        mpfr_clear(bb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

