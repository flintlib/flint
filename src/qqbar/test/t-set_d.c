/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "double_extras.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_set_d, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        double x;
        qqbar_t z;
        arb_t a, b;
        int ok;

        qqbar_init(z);
        arb_init(a);
        arb_init(b);

        x = d_randtest_special(state, -1100, 1100);

        ok =  qqbar_set_d(z, x);

        qqbar_get_arb(a, z, 53);

        arb_set_d(b, x);

        if (ok)
        {
            if (!arb_equal(a, b))
            {
                flint_printf("FAIL!\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            if (arb_is_finite(b))
            {
                flint_printf("FAIL!\n");
                flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
                flint_printf("a = "); arb_print(a); flint_printf("\n\n");
                flint_printf("b = "); arb_print(b); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        qqbar_clear(z);
    }

    TEST_FUNCTION_END(state);
}
