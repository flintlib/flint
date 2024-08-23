/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ctx_tau_set, state)
{
    slong iter;

    /* Test: various matrix equalities */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 10);
        arb_t pi;
        acb_mat_t tau;
        arb_mat_t N, id;
        acb_theta_ctx_tau_t ctx;

        arb_init(pi);
        acb_mat_init(tau, g, g);
        arb_mat_init(N, g, g);
        arb_mat_init(id, g, g);
        acb_theta_ctx_tau_init(ctx, g);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        acb_theta_ctx_tau_set(ctx, tau, prec);
        arb_mat_one(id);

        if (!acb_mat_equal(acb_theta_ctx_tau(ctx), tau))
        {
            flint_printf("FAIL (tau)\n");
            flint_abort();
        }

        arb_mat_mul(N, acb_theta_ctx_y(ctx), acb_theta_ctx_yinv(ctx), prec);
        if (!arb_mat_overlaps(N, id))
        {
            flint_printf("FAIL (Yinv)\n");
            arb_mat_printd(N, 5);
            flint_abort();
        }

        if (g > 1)
        {
            arb_mat_mul(N, acb_theta_ctx_cho(ctx), acb_theta_ctx_choinv(ctx), prec);
            if (!arb_mat_overlaps(N, id))
            {
                flint_printf("FAIL (choinv)\n");
                flint_abort();
            }

            arb_const_pi(pi, prec);
            arb_mat_transpose(N, acb_theta_ctx_cho(ctx));
            arb_mat_mul(N, N, acb_theta_ctx_cho(ctx), prec);
            arb_mat_scalar_div_arb(N, N, pi, prec);
            if (!arb_mat_overlaps(N, acb_theta_ctx_y(ctx)))
            {
                flint_printf("FAIL (cho)\n");
                flint_abort();
            }
        }

        arb_clear(pi);
        acb_mat_clear(tau);
        arb_mat_clear(N);
        arb_mat_clear(id);
        acb_theta_ctx_tau_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
