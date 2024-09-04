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

TEST_FUNCTION_START(acb_theta_sum_a0_tilde, state)
{
    slong iter;

    /* Test: matches sum_all_tilde */
    for (iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong mprec = 100 + n_randint(state, 100);
        slong prec = mprec + 50;
        acb_mat_t tau;
        acb_ptr z;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx;
        arb_ptr d;
        acb_ptr th1, th2;
        ulong a;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        acb_theta_ctx_z_init(ctx, g);
        d = _arb_vec_init(n);
        th1 = _acb_vec_init(n);
        th2 = _acb_vec_init(n * n);

        acb_siegel_randtest_compact(tau, state, 0, prec);
        acb_siegel_randtest_vec(z, state, g, prec);
        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx, z, ctx_tau, prec);
        acb_theta_agm_distances(d, z, 1, tau, prec);

        /* Call sum_a0_tilde at precision mprec, test against sum_all_tilde */
        acb_theta_sum_a0_tilde(th1, ctx, 1, ctx_tau, d, mprec);
        acb_theta_sum_all_tilde(th2, ctx, 1, ctx_tau, d, prec);
        for (a = 0; a < n; a++)
        {
            acb_set(&th2[a], &th2[a << g]);
        }

        if (!_acb_vec_overlaps(th1, th2, n))
        {
            flint_printf("FAIL\n");
            flint_printf("g=%wd\n", g);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("th1: ");
            _acb_vec_printd(th1, n, 5);
            flint_printf("th2: ");
            _acb_vec_printd(th2, n, 5);
            flint_printf("Difference: ");
            _acb_vec_sub(th1, th1, th2, n, prec);
            _acb_vec_printd(th1, n, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx);
        _arb_vec_clear(d, n);
        _acb_vec_clear(th1, n);
        _acb_vec_clear(th2, n * n);
    }

    TEST_FUNCTION_END(state);
}
