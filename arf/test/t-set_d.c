/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/double_extras.h"
#include "arf.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("set_d....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        double x;
        arf_t y, z;
        mpfr_t m;

        arf_init(y);
        arf_init(z);
        mpfr_init2(m, 53);

        x = d_randtest_special(state, -1200, 1200);
        arf_set_d(y, x);
        mpfr_set_d(m, x, MPFR_RNDN);
        arf_set_mpfr(z, m);

        if (!arf_equal(y, z) || !arf_equal_d(y, x))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = %.17g\n\n", x);
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(y);
        arf_clear(z);
        mpfr_clear(m);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

