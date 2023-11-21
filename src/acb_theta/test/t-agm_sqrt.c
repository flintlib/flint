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

TEST_FUNCTION_START(acb_theta_agm_sqrt, state)
{
    slong iter;

    /* Test:
       - if nonzero, value of square root should agree; precision remains high
       - no abort on wrong values of rt_low */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 1000);
        slong mag_bits = n_randint(state, 4);
        slong lowprec = 10 + n_randint(state, 10);
        slong delta = (1 << mag_bits) + 10;
        acb_t rt, x, rt_low, t;
        arf_t err;

        acb_init(rt);
        acb_init(x);
        acb_init(rt_low);
        acb_init(t);
        arf_init(err);

        acb_randtest_precise(rt, state, prec, mag_bits);
        if (iter % 10 == 0)
        {
            acb_union(rt, rt, x, prec);
        }

        acb_sqr(x, rt, prec);
        acb_one(t);
        acb_mul_2exp_si(t, t, -lowprec);
        acb_add_si(t, t, 1, lowprec);
        acb_mul(rt_low, rt, t, lowprec);

        acb_theta_agm_sqrt(t, x, rt_low, 1, prec);
        acb_get_rad_ubound_arf(err, t, prec);

        if (!acb_contains(t, rt))
        {
            flint_printf("FAIL (value)\n");
            acb_printd(rt, 5);
            flint_printf("\n");
            acb_printd(t, 5);
            flint_printf("\n");
            acb_printd(rt_low, 5);
            flint_printf("\n");
            flint_abort();
        }

        if (!acb_is_finite(t))
        {
            flint_printf("FAIL (infinite)\n");
            flint_abort();
        }

        if (!acb_contains_zero(rt) && (arf_cmp_2exp_si(err, -prec + delta) > 0))
        {
            flint_printf("FAIL (precision)\n");
            flint_printf("prec = %wd, result:\n", prec, mag_bits);
            acb_printd(t, 10);
            flint_printf("\nrt_low:\n");
            acb_printd(rt_low, 10);
            flint_printf("\n");
            flint_abort();
        }

        if (iter % 10 == 1)
        {
            acb_randtest(rt_low, state, prec, mag_bits);
            acb_theta_agm_sqrt(t, x, rt_low, 1, prec);
        }

        acb_clear(rt);
        acb_clear(x);
        acb_clear(rt_low);
        acb_clear(t);
        arf_clear(err);
    }

    TEST_FUNCTION_END(state);
}
