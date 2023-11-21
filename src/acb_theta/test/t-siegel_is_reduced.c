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

TEST_FUNCTION_START(acb_theta_siegel_is_reduced, state)
{
    slong iter;

    /* Test: correct values on some matrices */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = ACB_THETA_LOW_PREC;
        slong tol_exp = -10;
        slong j = n_randint(state, g);
        slong k = n_randint(state, g);
        acb_mat_t tau;

        acb_mat_init(tau, g, g);

        acb_mat_onei(tau);
        if (!acb_siegel_is_reduced(tau, tol_exp, prec))
        {
            flint_printf("FAIL (1)\n");
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_add_si(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k), 1, prec);
        acb_set(acb_mat_entry(tau, k, j), acb_mat_entry(tau, j, k));
        if (acb_siegel_is_reduced(tau, tol_exp, prec))
        {
            flint_printf("FAIL (2)\n");
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_mat_onei(tau);
        acb_mul_2exp_si(acb_mat_entry(tau, j, j), acb_mat_entry(tau, j, j), -1);
        if (acb_siegel_is_reduced(tau, tol_exp, prec))
        {
            flint_printf("FAIL (3)\n");
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
    }

    TEST_FUNCTION_END(state);
}
