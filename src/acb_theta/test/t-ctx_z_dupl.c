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

TEST_FUNCTION_START(acb_theta_ctx_z_dupl, state)
{
    slong iter;

    /* Test: matches with acb_theta_ctx_z_set with doubled input */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx1, ctx2;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        acb_theta_ctx_z_init(ctx1, g);
        acb_theta_ctx_z_init(ctx2, g);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx1, z, ctx_tau, prec);
        acb_theta_ctx_z_dupl(ctx1, prec);

        acb_theta_ctx_tau_dupl(ctx_tau, prec);
        _acb_vec_scalar_mul_2exp_si(z, z, g, 1);
        acb_theta_ctx_z_set(ctx2, z, ctx_tau, prec);

        if (!acb_theta_ctx_z_overlaps(ctx1, ctx2))
        {
            flint_printf("FAIL\n");
            _acb_vec_printd(acb_theta_ctx_z(ctx1), g, 5);
            _acb_vec_printd(acb_theta_ctx_z(ctx2), g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_z(ctx1), g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_z(ctx2), g, 5);
            acb_printd(acb_theta_ctx_c(ctx1), 5);
            flint_printf("\n");
            acb_printd(acb_theta_ctx_c(ctx2), 5);
            flint_printf("\n");
            _arb_vec_printd(acb_theta_ctx_r(ctx1), g, 5);
            _arb_vec_printd(acb_theta_ctx_r(ctx2), g, 5);
            if (g > 1)
            {
                _acb_vec_printd(acb_theta_ctx_exp_2z(ctx1), g, 5);
                _acb_vec_printd(acb_theta_ctx_exp_2z(ctx2), g, 5);
                _acb_vec_printd(acb_theta_ctx_exp_z_inv(ctx1), g, 5);
                _acb_vec_printd(acb_theta_ctx_exp_z_inv(ctx2), g, 5);
                _acb_vec_printd(acb_theta_ctx_exp_2z_inv(ctx1), g, 5);
                _acb_vec_printd(acb_theta_ctx_exp_2z_inv(ctx2), g, 5);
                _arb_vec_printd(acb_theta_ctx_v(ctx1), g, 5);
                _arb_vec_printd(acb_theta_ctx_v(ctx2), g, 5);
            }
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx1);
        acb_theta_ctx_z_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
