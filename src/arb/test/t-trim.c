/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_trim, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y;
        slong acc1, acc2;
        int accuracy_ok;

        arb_init(x);
        arb_init(y);

        arb_randtest_special(x, state, 1000, 100);

        arb_trim(y, x);

        if (!arb_contains(y, x))
        {
            flint_printf("FAIL (containment):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(x);
        acc2 = arb_rel_accuracy_bits(y);

        accuracy_ok = (acc1 < 0 && acc2 <= acc1) || (acc2 >= acc1 - 1);

        if (!accuracy_ok)
        {
            flint_printf("FAIL (accuracy):\n\n");
            flint_printf("x: %wd, y = %wd\n\n", acc1, acc2);
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_trim(x, x);

        if (!arb_equal(y, x))
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
