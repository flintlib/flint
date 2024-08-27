/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ctx_tau_copy, state)
{
    slong iter;

    /* Test: copy overlaps original */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 10);
        acb_mat_t tau;
        acb_theta_ctx_tau_t ctx, copy;

        acb_mat_init(tau, g, g);
        acb_theta_ctx_tau_init(ctx, g);
        acb_theta_ctx_tau_init(copy, g);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_ctx_tau_set(ctx, tau, prec);
        acb_theta_ctx_tau_copy(copy, ctx);

        if (!acb_theta_ctx_tau_overlaps(ctx, copy))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_theta_ctx_tau_clear(ctx);
        acb_theta_ctx_tau_clear(copy);
    }

    TEST_FUNCTION_END(state);
}
