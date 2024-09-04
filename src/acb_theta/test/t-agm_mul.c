/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_agm_mul, state)
{
    slong iter;

    /* Test: duplication formula using sum_a0_tilde */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong mprec = 100 + n_randint(state, 200);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr z;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx0, ctx;
        arb_ptr d0, d;
        acb_ptr th2, th_dupl, test;
        slong j;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        acb_theta_ctx_z_init(ctx0, g);
        acb_theta_ctx_z_init(ctx, g);
        th2 = _acb_vec_init(2 * n);
        th_dupl = _acb_vec_init(2 * n);
        test = _acb_vec_init(2 * n);
        d0 = _arb_vec_init(n);
        d = _arb_vec_init(n);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(z + g, state, tau, 0, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx0, z, ctx_tau, prec);
        acb_theta_ctx_z_set(ctx, z + g, ctx_tau, prec);
        acb_theta_dist_a0(d0, z, tau, prec);
        acb_theta_dist_a0(d, z + g, tau, prec);

        /* Make test vector using sum_a0_tilde squared */
        acb_theta_sum_a0_tilde(test, ctx0, 1, ctx_tau, d0, prec);
        acb_theta_sum_a0_tilde(test + n, ctx, 1, ctx_tau, d, prec);
        for (j = 0; j < 2 * n; j++)
        {
            acb_sqr(&test[j], &test[j], prec);
        }

        /* Duplicate to get input of agm_mul */
        acb_theta_ctx_tau_dupl(ctx_tau, prec);
        acb_theta_ctx_z_dupl(ctx0, prec);
        acb_theta_ctx_z_dupl(ctx, prec);
        _arb_vec_scalar_mul_2exp_si(d0, d0, n, 1);
        _arb_vec_scalar_mul_2exp_si(d, d, n, 1);
        acb_theta_sum_a0_tilde(th_dupl, ctx0, 1, ctx_tau, d0, prec);
        acb_theta_sum_a0_tilde(th_dupl + n, ctx, 1, ctx_tau, d, prec);

        /* Call agm_mul at precision mprec and compare with test */
        acb_theta_agm_mul(th2, th_dupl, th_dupl, g, mprec);
        acb_theta_agm_mul(th2 + n, th_dupl, th_dupl + n, g, mprec);
        _acb_vec_scalar_mul_2exp_si(th2, th2, 2 * n, g);

        if (!_acb_vec_overlaps(test, th2, 2 * n))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("g = %wd, prec = %wd, tau, z:\n", g, prec);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, 2 * g, 5);
            flint_printf("th2:\n");
            _acb_vec_printd(th2, 2 * n, 5);
            flint_printf("th_dupl:\n");
            _acb_vec_printd(th_dupl, 2 * n, 5);
            flint_printf("test:\n");
            _acb_vec_printd(test, 2 * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx0);
        acb_theta_ctx_z_clear(ctx);
        _acb_vec_clear(th2, 2 * n);
        _acb_vec_clear(th_dupl, 2 * n);
        _acb_vec_clear(test, 2 * n);
        _arb_vec_clear(d0, n);
        _arb_vec_clear(d, n);
    }

    TEST_FUNCTION_END(state);
}
