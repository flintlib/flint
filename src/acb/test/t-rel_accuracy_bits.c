/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_rel_accuracy_bits, state)
{
    slong iter;

    /* test aliasing of c and a */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        acb_t z;
        slong a1, a2;

        arb_init(x);
        acb_init(z);

        arb_randtest_special(x, state, 1 + n_randint(state, 200), 1 + n_randint(state, 200));
        acb_set_arb(z, x);

        a1 = arb_rel_accuracy_bits(x);
        a2 = acb_rel_accuracy_bits(z);

        if (a1 != a2)
        {
            flint_printf("FAIL: acb != arb\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_printf("a1 = %wd, a2 = %wd\n\n", a1, a2);
            flint_abort();
        }

        acb_randtest_special(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 200));

        a1 = acb_rel_accuracy_bits(z);

        if (n_randint(state, 2))
            arf_swap(arb_midref(acb_realref(z)), arb_midref(acb_imagref(z)));

        if (n_randint(state, 2))
            mag_swap(arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));

        a2 = acb_rel_accuracy_bits(z);

        if (a1 != a2)
        {
            flint_printf("FAIL: swapping\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_printf("a1 = %wd, a2 = %wd\n\n", a1, a2);
            flint_abort();
        }

        acb_randtest_special(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 200));

        if (arf_cmpabs(arb_midref(acb_realref(z)), arb_midref(acb_imagref(z))) >= 0)
            arf_set(arb_midref(x), arb_midref(acb_realref(z)));
        else
            arf_set(arb_midref(x), arb_midref(acb_imagref(z)));

        if (mag_cmp(arb_radref(acb_realref(z)), arb_radref(acb_imagref(z))) >= 0)
            mag_set(arb_radref(x), arb_radref(acb_realref(z)));
        else
            mag_set(arb_radref(x), arb_radref(acb_imagref(z)));

        a1 = acb_rel_accuracy_bits(z);
        a2 = arb_rel_accuracy_bits(x);

        if (a1 != a2)
        {
            flint_printf("FAIL: acb != arb (2)\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_printf("a1 = %wd, a2 = %wd\n\n", a1, a2);
            flint_abort();
        }

        arb_clear(x);
        acb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
