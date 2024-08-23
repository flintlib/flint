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

TEST_FUNCTION_START(acb_theta_ctx_exp_inv, state)
{
    slong iter;

    /* Test: correct value, never indeterminate */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_t z, exp, exp_inv, test;
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 4);
        int is_real = iter % 2;

        acb_init(z);
        acb_init(exp);
        acb_init(exp_inv);
        acb_init(test);

        acb_randtest(z, state, prec, mag_bits);
        if (is_real)
        {
            arb_zero(acb_imagref(z));
        }
        acb_exp_pi_i(exp, z, prec);
        acb_inv(exp_inv, exp, prec);
        acb_theta_ctx_exp_inv(test, exp, z, is_real, prec);

        if (!acb_overlaps(test, exp_inv)
            || !acb_is_finite(test))
        {
            flint_printf("FAIL\n");
            acb_printd(test, 5);
            flint_printf("\n");
            acb_printd(exp_inv, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mul(test, test, exp, prec);
        acb_sub_si(test, test, 1, prec);

        if (!acb_contains_zero(test))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(exp);
        acb_clear(exp_inv);
        acb_clear(test);
    }

    TEST_FUNCTION_END(state);
}
