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

TEST_FUNCTION_START(arb_hypgeom_gamma_taylor, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, s1, s2, a, b;
        slong prec, ebits, prec2;
        int success, success2, alias, reciprocal;

        if (n_randint(state, 10) == 0)
            prec = 2 + n_randint(state, 4000);
        else
            prec = 2 + n_randint(state, 300);

        if (n_randint(state, 10) == 0)
            ebits = 100;
        else
            ebits = 10;

        prec2 = prec + 1 + n_randint(state, 30);

        arb_init(x);
        arb_init(s1);
        arb_init(s2);
        arb_init(a);
        arb_init(b);

        arb_randtest(x, state, prec, ebits);
        arb_randtest(s1, state, prec, 10);
        arb_randtest(s2, state, prec, 10);
        alias = n_randint(state, 2);
        reciprocal = n_randint(state, 2);

        if (alias)
        {
            success = arb_hypgeom_gamma_taylor(s1, x, reciprocal, prec);
        }
        else
        {
            arb_set(s1, x);
            success = arb_hypgeom_gamma_taylor(s1, s1, reciprocal, prec);
        }

        if (success)
        {
            /* printf("%ld\n", iter); */

            /* Compare with Stirling series algorithm. */
            arb_hypgeom_gamma_stirling(s2, x, reciprocal, prec);

            if (!arb_overlaps(s1, s2))
            {
                flint_printf("FAIL\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); arb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
                arb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }

            /* Compare with different level of precision. */
            success2 = arb_hypgeom_gamma_taylor(s2, x, reciprocal, prec2);

            if (success2 && !arb_overlaps(s1, s2))
            {
                flint_printf("FAIL (2)\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); arb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
                arb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }

            arf_set_mag(arb_midref(a), arb_radref(x));
            arf_set_mag(arb_midref(b), arb_radref(x));

            arb_sub_arf(a, a, arb_midref(x), prec + 30);
            arb_neg(a, a);

            arb_add_arf(b, b, arb_midref(x), prec + 30);

            success2 = arb_hypgeom_gamma_taylor(s2, a, reciprocal, prec2);

            if (success2 && !arb_overlaps(s1, s2))
            {
                flint_printf("FAIL (3)\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); arb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
                arb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }

            success2 = arb_hypgeom_gamma_taylor(s2, b, reciprocal, prec2);

            if (success2 && !arb_overlaps(s1, s2))
            {
                flint_printf("FAIL (4)\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); arb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
                arb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }

            arb_add(a, a, b, prec + 30);
            arb_mul_2exp_si(a, a, -1);

            success2 = arb_hypgeom_gamma_taylor(s2, b, reciprocal, prec2);

            if (success2 && !arb_overlaps(s1, s2))
            {
                flint_printf("FAIL (5)\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); arb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); arb_printn(s2, 1000, 0); flint_printf("\n\n");
                arb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); arb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(x);
        arb_clear(s1);
        arb_clear(s2);
        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
