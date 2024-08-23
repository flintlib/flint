/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ctx_sqr_inv, state)
{
    slong iter;

    /* Test: correct value, never indeterminate */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_t z, exp, sqr, inv, sqr_inv, test;
        slong prec = 10 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 4);
        int is_real = iter % 2;

        acb_init(z);
        acb_init(exp);
        acb_init(sqr);
        acb_init(inv);
        acb_init(sqr_inv);
        acb_init(test);

        acb_randtest(z, state, prec, mag_bits);
        if (is_real)
        {
            arb_zero(acb_imagref(z));
        }
        acb_exp_pi_i(exp, z, prec);
        acb_sqr(sqr, exp, prec);

        acb_neg(z, z);
        acb_exp_pi_i(inv, z, prec);
        acb_mul_2exp_si(z, z, 1);
        acb_exp_pi_i(sqr_inv, z, prec);

        acb_theta_ctx_sqr_inv(test, inv, sqr, is_real, prec);

        if (!acb_overlaps(test, sqr_inv)
            || !acb_is_finite(test))
        {
            flint_printf("FAIL\n");
            acb_printd(test, 5);
            flint_printf("\n");
            acb_printd(sqr_inv, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mul(test, test, sqr, prec);
        acb_sub_si(test, test, 1, prec);

        if (!acb_contains_zero(test))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(exp);
        acb_clear(sqr);
        acb_clear(inv);
        acb_clear(sqr_inv);
        acb_clear(test);
    }

    TEST_FUNCTION_END(state);
}
