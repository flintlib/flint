/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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
        acb_ptr th;
        slong k;
        acb_t r, s;

        th = _acb_vec_init(n2);
        acb_init(r);
        acb_init(s);

        for (k = 0; k < n2; k++)
        {
            acb_randtest_precise(&th[k], state, prec, mag_bits);
        }

        acb_theta_g2_chi5(r, th, prec);
        acb_sqr(r, r, prec);
        _acb_vec_sqr(th, th, n2, prec);
        acb_theta_g2_chi10(s, th, prec);

        if (!acb_overlaps(r, s))
        {
            flint_printf("FAIL\n");
            acb_printd(r, 10);
            flint_printf("\n");
            acb_printd(s, 10);
            flint_printf("\n");
            flint_abort();
        }

        _acb_vec_clear(th, n2);
        acb_clear(r);
        acb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
