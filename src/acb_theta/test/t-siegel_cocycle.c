/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_siegel_cocycle, state)
{
    slong iter;

    /* Test: chain rule */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 10);
        fmpz_mat_t m1, m2, m3;
        acb_mat_t tau1, tau2;
        acb_mat_t c1, c2, c3, t;

        fmpz_mat_init(m1, 2 * g, 2 * g);
        fmpz_mat_init(m2, 2 * g, 2 * g);
        fmpz_mat_init(m3, 2 * g, 2 * g);
        acb_mat_init(tau1, g, g);
        acb_mat_init(tau2, g, g);
        acb_mat_init(c1, g, g);
        acb_mat_init(c2, g, g);
        acb_mat_init(c3, g, g);
        acb_mat_init(t, g, g);

        acb_siegel_randtest(tau1, state, prec, mag_bits);
        sp2gz_randtest(m1, state, mag_bits);
        sp2gz_randtest(m2, state, mag_bits);

        /* Test: chain rule */
        acb_siegel_cocycle(c1, m1, tau1, prec);
        acb_siegel_transform(tau2, m1, tau1, prec);
        acb_siegel_cocycle(c2, m2, tau2, prec);
        fmpz_mat_mul(m3, m2, m1);
        acb_siegel_cocycle(c3, m3, tau1, prec);
        acb_mat_mul(t, c2, c1, prec);

        if (!acb_mat_overlaps(t, c3))
        {
            flint_printf("FAIL\n\n");
            acb_mat_printd(c3, 10);
            flint_printf("\n");
            acb_mat_printd(t, 10);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(m1);
        fmpz_mat_clear(m2);
        fmpz_mat_clear(m3);
        acb_mat_clear(tau1);
        acb_mat_clear(tau2);
        acb_mat_clear(c1);
        acb_mat_clear(c2);
        acb_mat_clear(c3);
        acb_mat_clear(t);
    }

    TEST_FUNCTION_END(state);
}
