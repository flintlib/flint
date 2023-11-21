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
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_2f1_continuation, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z1, z2, f1, f2, g1, g2, h1, h2, aa, bb, cc;
        mag_t d0, d1, dt;

        slong prec;
        int regularized, ebits;

        acb_init(a); acb_init(b); acb_init(c);
        acb_init(aa); acb_init(bb); acb_init(cc);
        acb_init(z1); acb_init(z2);
        acb_init(f1); acb_init(f2);
        acb_init(g1); acb_init(g2);
        acb_init(h1); acb_init(h2);
        mag_init(d0); mag_init(d1); mag_init(dt);

        prec = 2 + n_randint(state, 300);
        ebits = 10;
        regularized = n_randint(state, 2);

        acb_randtest_param(a, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits / 2));
        acb_randtest_param(b, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits / 2));
        acb_randtest_param(c, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits / 2));
        acb_randtest(h1, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        acb_randtest(h2, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));

        do {
            int left, upper, lower;

            acb_randtest_param(z1, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
            acb_randtest_param(z2, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));

            /* we test both convergent and non-convergent cases, but
               try to be more efficient by generating more convergent cases */
            if (n_randint(state, 2))
            {
                acb_sub_ui(aa, z1, 1, prec);
                acb_get_mag(d0, z1);
                acb_get_mag(d1, aa);
                acb_get_mag(dt, z2);

                if (mag_cmp(dt, d0) >= 0 || mag_cmp(dt, d1) >= 0)
                    continue;
            }

            acb_add(z2, z1, z2, prec);

            /* for the test, don't cross the branch cut */
            acb_sub_ui(aa, z1, 1, prec);
            acb_sub_ui(bb, z2, 1, prec);

            left = arb_is_negative(acb_realref(aa)) && arb_is_negative(acb_realref(bb));
            upper = arb_is_positive(acb_imagref(aa)) && arb_is_positive(acb_imagref(bb));
            lower = arb_is_nonpositive(acb_imagref(aa)) && arb_is_nonpositive(acb_imagref(bb));

            if (left || upper || lower)
                break;
        } while (1);

        acb_add_ui(aa, a, 1, prec);
        acb_add_ui(bb, b, 1, prec);
        acb_add_ui(cc, c, 1, prec);

        acb_hypgeom_2f1(f1, a, b, c, z1, regularized, prec);
        acb_hypgeom_2f1(f2, aa, bb, cc, z1, regularized, prec);
        acb_mul(f2, f2, a, prec);
        acb_mul(f2, f2, b, prec);
        if (!regularized)
            acb_div(f2, f2, c, prec);

        acb_hypgeom_2f1_continuation(h1, h2, a, b, c, z1, z2, f1, f2, prec);

        if (acb_is_finite(h1) || acb_is_finite(h2))
        {
            acb_hypgeom_2f1(g1, a, b, c, z2, regularized, prec);
            acb_hypgeom_2f1(g2, aa, bb, cc, z2, regularized, prec);
            acb_mul(g2, g2, a, prec);
            acb_mul(g2, g2, b, prec);
            if (!regularized)
                acb_div(g2, g2, c, prec);

            if (!acb_overlaps(g1, h1) || !acb_overlaps(g2, h2))
            {
                flint_printf("FAIL: consistency\n\n");
                flint_printf("regularized = %d, prec = %wd\n\n", regularized, prec);
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n\n");
                flint_printf("F(a,b,c,z1) and F'(a,b,c,z1):\n");
                flint_printf("f1 = "); acb_printd(f1, 30); flint_printf("\n\n");
                flint_printf("f2 = "); acb_printd(f2, 30); flint_printf("\n\n");
                flint_printf("F(a,b,c,z2) and F'(a,b,c,z2):\n");
                flint_printf("g1 = "); acb_printd(g1, 30); flint_printf("\n\n");
                flint_printf("g2 = "); acb_printd(g2, 30); flint_printf("\n\n");
                flint_printf("Computed F and F':\n");
                flint_printf("h1 = "); acb_printd(h1, 30); flint_printf("\n\n");
                flint_printf("h2 = "); acb_printd(h2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c);
        acb_clear(aa); acb_clear(bb); acb_clear(cc);
        acb_clear(z1); acb_clear(z2);
        acb_clear(f1); acb_clear(f2);
        acb_clear(g1); acb_clear(g2);
        acb_clear(h1); acb_clear(h2);
        mag_clear(d0); mag_clear(d1); mag_clear(dt);
    }

    TEST_FUNCTION_END(state);
}
