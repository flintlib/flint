/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_chi5, state)
{
    slong iter;

    /* Test: square is chi10 */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n2 = 1 << (2 * g);
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 10);
        acb_ptr th, mf;
        slong k;
        acb_t r;

        th = _acb_vec_init(n2);
        mf = _acb_vec_init(4);
        acb_init(r);

        for (k = 0; k < n2; k++)
        {
            acb_randtest_precise(&th[k], state, prec, mag_bits);
        }

        acb_theta_g2_chi5(r, th, prec);
        acb_sqr(r, r, prec);
        _acb_vec_sqr(th, th, n2, prec);
        acb_theta_g2_even_weight(&mf[0], &mf[1], &mf[2], &mf[3], th, prec);

        if (!acb_overlaps(r, &mf[2])
            || !_acb_vec_is_finite(mf, 4)
            || !acb_is_finite(r))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(&mf[2], 10);
            flint_printf("\n");
            flint_abort();
        }

        _acb_vec_clear(th, n2);
        _acb_vec_clear(mf, 4);
        acb_clear(r);
    }

    TEST_FUNCTION_END(state);
}
