/*
    Copyright (C) 2023 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_minmax, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z1, z2, w1, w2;
        slong prec;

        arb_init(x);
        arb_init(y);
        arb_init(z1);
	arb_init(z2);
	arb_init(w1);
	arb_init(w2);

	prec = 2 + n_randint(state, 200);

	arb_randtest_special(x, state, 1 + n_randint(state, 200), 10);
        arb_randtest_special(y, state, 1 + n_randint(state, 200), 10);

	arb_min(w1, x, y, prec);
	arb_max(w2, x, y, prec);
        arb_minmax(z1, z2, x, y, prec);

        if (!arb_equal(w1, z1) || !arb_equal(w2, z2))
        {
            flint_printf("FAIL: same as min and max\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z1 = "); arb_print(z1); flint_printf("\n\n");
	    flint_printf("z2 = "); arb_print(z2); flint_printf("\n\n");
	    flint_printf("w1 = "); arb_print(w1); flint_printf("\n\n");
	    flint_printf("w2 = "); arb_print(w2); flint_printf("\n\n");
            flint_abort();
        }

        /* aliasing */
        {
            int alias;

            if (n_randint(state, 2))
            {
		arb_minmax(x, y, x, y, prec);
                alias = arb_equal(x, z1) && arb_equal(y, z2);
            }
            else
            {
                arb_minmax(y, x, x, y, prec);
                alias = arb_equal(y, z1) && arb_equal(x, z2);
            }

            if (!alias)
            {
                flint_printf("FAIL: aliasing\n\n");
                flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                flint_printf("y = "); arb_print(y); flint_printf("\n\n");
                flint_printf("z1 = "); arb_print(z1); flint_printf("\n\n");
		flint_printf("z2 = "); arb_print(z2); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z1);
	arb_clear(z2);
	arb_clear(w1);
	arb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}
