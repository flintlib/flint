/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_gamma_taylor, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, s1, s2, a, b;
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

        acb_init(x);
        acb_init(s1);
        acb_init(s2);
        acb_init(a);
        acb_init(b);

        acb_randtest(x, state, prec, ebits);
        acb_randtest(s1, state, prec, 10);
        acb_randtest(s2, state, prec, 10);
        alias = n_randint(state, 2);
        reciprocal = n_randint(state, 2);

        if (alias)
        {
            success = acb_hypgeom_gamma_taylor(s1, x, reciprocal, prec);
        }
        else
        {
            acb_set(s1, x);
            success = acb_hypgeom_gamma_taylor(s1, s1, reciprocal, prec);
        }

        if (success)
        {
            /* printf("%ld\n", iter); */

            /* Compare with Stirling series algorithm. */
            acb_hypgeom_gamma_stirling(s2, x, reciprocal, prec);

            if (!acb_overlaps(s1, s2))
            {
                flint_printf("FAIL\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); acb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); acb_printn(s2, 1000, 0); flint_printf("\n\n");
                acb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); acb_printd(s1, 1000); flint_printf("\n\n");
                flint_abort();
            }

            /* Compare with different level of precision. */
            success2 = acb_hypgeom_gamma_taylor(s2, x, reciprocal, prec2);

            if (success2 && !acb_overlaps(s1, s2))
            {
                flint_printf("FAIL (2)\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); acb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); acb_printn(s2, 1000, 0); flint_printf("\n\n");
                acb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }

            acb_get_mid(a, x);

            if (n_randint(state, 2))
            {
                arf_set_mag(arb_midref(acb_realref(b)), arb_radref(acb_realref(x)));
                arf_set_mag(arb_midref(acb_imagref(b)), arb_radref(acb_imagref(x)));

                if (n_randint(state, 2))
                    acb_neg(b, b);
                if (n_randint(state, 2))
                    acb_conj(b, b);

                acb_add(a, a, b, prec2);
            }

            success2 = acb_hypgeom_gamma_taylor(s2, a, reciprocal, prec2);

            if (success2 && !acb_overlaps(s1, s2))
            {
                flint_printf("FAIL (3)\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("x = "); acb_printn(x, 1000, 0); flint_printf("\n\n");
                flint_printf("s1 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_printf("s2 = "); acb_printn(s2, 1000, 0); flint_printf("\n\n");
                acb_sub(s1, s1, s2, prec2);
                flint_printf("s1 - s2 = "); acb_printn(s1, 1000, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(x);
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(a);
        acb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
