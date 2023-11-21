/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_transform_sqrtdet, state)
{
    slong iter;

    /* Test: square of sqrtdet is det */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        acb_mat_t tau;
        acb_t r, t;
        slong prec = 2 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 4);

        acb_mat_init(tau, g, g);
        acb_init(r);
        acb_init(t);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_transform_sqrtdet(r, tau, prec);
        acb_sqr(r, r, prec);
        acb_mat_det(t, tau, prec);

        if (!acb_overlaps(r, t))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(t, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_clear(r);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
