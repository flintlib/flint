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

TEST_FUNCTION_START(acb_theta_ctx_set_t, state)
{
    slong iter;

    /* Test: corresponds to acb_theta_ctx_set_z for 0, t, 2t, ... */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z, t, w;
        acb_theta_ctx_t ctx1, ctx2;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        w = _acb_vec_init(g);
        acb_theta_ctx_init(ctx1, 6, g);
        acb_theta_ctx_init(ctx2, 6, g);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g, prec);
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec);
        }

        /* Sometimes force z to be real or zero, or t to be zero */
        if (iter % 3 == 0)
        {
            for (k = 0; k < g; k++)
            {
                arb_zero(acb_imagref(&z[k]));
            }
        }
        if (iter % 4 == 0)
        {
            _acb_vec_zero(t, g);
        }
        if (iter % 6 == 0)
        {
            _acb_vec_zero(z, g);
        }

        acb_theta_ctx_set_tau(ctx1, tau, prec);
        acb_theta_ctx_set_z_ql(ctx1, z, prec);
        acb_theta_ctx_set_t(ctx1, t, prec);

        acb_theta_ctx_set_tau(ctx2, tau, prec);
        acb_theta_ctx_set_z(ctx2, w, 0, prec);
        if (!_acb_vec_is_zero(z, g))
        {
            acb_theta_ctx_set_z(ctx2, z, 3, prec);
        }
        if (!_acb_vec_is_zero(t, g))
        {
            acb_theta_ctx_set_z(ctx2, t, 1, prec);
            _acb_vec_scalar_mul_2exp_si(w, t, g, 1);
            acb_theta_ctx_set_z(ctx2, w, 2, prec);
        }
        if (!_acb_vec_is_zero(z, g) && !_acb_vec_is_zero(t, g))
        {
            _acb_vec_add(w, z, t, g, prec);
            acb_theta_ctx_set_z(ctx2, w, 4, prec);
            _acb_vec_add(w, w, t, g, prec);
            acb_theta_ctx_set_z(ctx2, w, 5, prec);
        }

        if (!_acb_vec_overlaps(acb_theta_ctx_exp_zs(ctx1), acb_theta_ctx_exp_zs(ctx2), 6 * g))
        {
            flint_printf("FAIL (exp_zs)\n");
            _acb_vec_printd(z, g, 5);
            _acb_vec_printd(t, g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx1), 6 * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx2), 6 * g, 5);
            flint_abort();
        }
        if (!_acb_vec_overlaps(acb_theta_ctx_exp_zs_inv(ctx1), acb_theta_ctx_exp_zs_inv(ctx2), 6 * g))
        {
            flint_printf("FAIL (exp_zs_inv)\n");
            _acb_vec_printd(z, g, 5);
            _acb_vec_printd(t, g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs_inv(ctx1), 6 * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs_inv(ctx2), 6 * g, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(w, g);
        acb_theta_ctx_clear(ctx1);
        acb_theta_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
