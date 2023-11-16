/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_lgamma, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t z, s1, s2, a, b;
        slong prec, ebits, prec2;

        prec = 2 + n_randint(state, 200);

        if (n_randint(state, 10) == 0)
            prec = 2 + n_randint(state, 1000);

        if (n_randint(state, 10) == 0)
            ebits = 100;
        else
            ebits = 10;
        ebits = 2;

        prec2 = 2 + n_randint(state, 200);

        arb_init(z);
        arb_init(s1);
        arb_init(s2);
        arb_init(a);
        arb_init(b);

        arb_randtest(z, state, prec, ebits);
        arb_randtest(s1, state, prec, 10);
        arb_randtest(s2, state, prec, 10);

        if (n_randint(state, 2))
        {
            arb_hypgeom_lgamma(s1, z, prec);
        }
        else
        {
            arb_set(s1, z);
            arb_hypgeom_lgamma(s1, s1, prec);
        }

        arb_add_ui(s2, z, 1, prec2);
        arb_hypgeom_lgamma(s2, s2, prec2);
        arb_log(a, z, prec2);
        arb_sub(s2, s2, a, prec2);

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("z = "); arb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
            arb_sub(s1, s1, s2, prec2);
            flint_printf("s1 - s2 = "); arb_printd(s1, 1000); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(a, z);
        mag_zero(arb_radref(a));

        if (n_randint(state, 2))
        {
            arf_set_mag(arb_midref(b), arb_radref(z));

            if (n_randint(state, 2))
                arb_neg(b, b);

            arb_add(a, a, b, prec);
        }

        arb_hypgeom_lgamma(s2, a, prec);

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL (2)\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("z = "); arb_printn(z, 1000, 0); flint_printf("\n\n");
            flint_printf("a = "); arb_printn(a, 1000, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
            arb_sub(s1, s1, s2, prec2);
            flint_printf("s1 - s2 = "); arb_printd(s1, 1000); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(z);
        arb_clear(s1);
        arb_clear(s2);
        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
