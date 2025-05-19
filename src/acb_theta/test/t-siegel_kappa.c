/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_siegel_kappa, state)
{
    slong iter;

    /* Test: compatibility with matrix multiplication, and compatible with sqr */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong bits = n_randint(state, 4);
        slong prec = 200;
        fmpz_mat_t m1, m2, m;
        acb_mat_t tau;
        acb_t s1, s2, s, t;
        slong kappa1, kappa2, kappa, e1, e2, e;
        ulong c1, c2, c;

        fmpz_mat_init(m1, 2 * g, 2 * g);
        fmpz_mat_init(m2, 2 * g, 2 * g);
        fmpz_mat_init(m, 2 * g, 2 * g);
        acb_mat_init(tau, g, g);
        acb_init(s1);
        acb_init(s2);
        acb_init(s);
        acb_init(t);

        sp2gz_randtest(m1, state, bits);
        sp2gz_randtest(m2, state, bits);
        acb_siegel_randtest_reduced(tau, state, prec, bits);

        fmpz_mat_mul(m, m2, m1);
        kappa = acb_siegel_kappa(s, m, tau, 0, prec);

        kappa1 = acb_siegel_kappa(s1, m1, tau, 0, prec);
        acb_siegel_transform(tau, m1, tau, prec);
        kappa2 = acb_siegel_kappa(s2, m2, tau, 0, prec);

        acb_theta_char_table(&c2, &e2, m2, 0, 0);
        acb_theta_char_table(&c1, &e1, m1, c2, 0);
        acb_theta_char_table(&c, &e, m, 0, 0);

        if (c1 != c
            || ((kappa1 + e1 + kappa2 + e2) % 4 != (kappa + e) % 4))
        {
            flint_printf("FAIL (characteristics)\n");
            flint_printf("c1 = %wd, c2 = %wd, c = %wd, e1 = %wd, e2 = %wd, e = %wd, kappa1 = %wd, kappa2 = %wd, kappa = %wd\n",
                c1, c2, c, e1, e2, e, kappa1, kappa2, kappa);
            flint_abort();
        }

        acb_mul(t, s1, s2, prec);
        if ((kappa1 + e1 + kappa2 + e2) % 8 != (kappa + e) % 8)
        {
            acb_neg(t, t);
        }
        if (!acb_overlaps(t, s)
            || !acb_is_finite(t)
            || !acb_is_finite(s))
        {
            flint_printf("FAIL (square roots)\n");
            flint_printf("s, t:\n");
            acb_printd(t, 5);
            flint_printf("\n");
            acb_printd(s, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_sqr(t, s2, prec);
        kappa = acb_siegel_kappa(s2, m2, tau, 1, prec);
        if ((kappa != kappa2 % 4)
            || !acb_overlaps(t, s2)
            || !acb_is_finite(t)
            || !acb_is_finite(s2))
        {
            flint_printf("FAIL (sqr)\n");
            flint_printf("kappa2 = %wd, kappa = %wd, s2, t:\n", kappa2, kappa);
            acb_printd(s2, 5);
            flint_printf("\n");
            acb_printd(t, 5);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(m1);
        fmpz_mat_clear(m2);
        fmpz_mat_clear(m);
        acb_mat_clear(tau);
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(s);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
