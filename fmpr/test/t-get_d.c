/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_d....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip */
    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, z;
        double y;

        fmpr_init(x);
        fmpr_init(z);

        fmpr_randtest_special(x, state, 53, 8);
        y = fmpr_get_d(x, FMPR_RND_DOWN);
        fmpr_set_d(z, y);

        if (!fmpr_equal(x, z))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = %.17g\n\n", y);
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

