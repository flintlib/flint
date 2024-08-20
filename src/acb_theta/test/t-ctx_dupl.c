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

TEST_FUNCTION_START(acb_theta_ctx_dupl, state)
{
    slong iter;

    /* Test: matches with acb_theta_ctx_set with doubled input */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z, t;
        acb_theta_ctx_t ctx1, ctx2;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
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
        acb_theta_ctx_dupl(ctx1, prec);

        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        _acb_vec_scalar_mul_2exp_si(z, z, g, 1);
        _acb_vec_scalar_mul_2exp_si(t, t, g, 1);
        acb_theta_ctx_set_tau(ctx2, tau, prec);
        acb_theta_ctx_set_z_ql(ctx2, z, prec);
        acb_theta_ctx_set_t(ctx2, t, prec);

        if (!acb_mat_overlaps(acb_theta_ctx_tau(ctx1), acb_theta_ctx_tau(ctx2)))
        {
            flint_printf("FAIL (tau)\n");
            acb_mat_printd(acb_theta_ctx_tau(ctx1), 5);
            acb_mat_printd(acb_theta_ctx_tau(ctx2), 5);
            flint_abort();
        }
        if (!arb_mat_overlaps(acb_theta_ctx_y(ctx1), acb_theta_ctx_y(ctx2)))
        {
            flint_printf("FAIL (Y)\n");
            flint_abort();
        }
        if (!arb_mat_overlaps(acb_theta_ctx_yinv(ctx1), acb_theta_ctx_yinv(ctx2)))
        {
            flint_printf("FAIL (Yinv)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_4(ctx1), acb_theta_ctx_exp_tau_div_4(ctx2)))
        {
            flint_printf("FAIL (exp_tau_div_4)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_2(ctx1), acb_theta_ctx_exp_tau_div_2(ctx2)))
        {
            flint_printf("FAIL (exp_tau_div_2)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau(ctx1), acb_theta_ctx_exp_tau(ctx2)))
        {
            flint_printf("FAIL (exp_tau)\n");
            flint_abort();
        }
        if (!_acb_vec_overlaps(acb_theta_ctx_exp_zs(ctx1), acb_theta_ctx_exp_zs(ctx2), 6 * g))
        {
            flint_printf("FAIL (exp_zs)\n");
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx1), 6 * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_zs(ctx2), 6 * g, 5);
            flint_abort();
        }
        if (!_acb_vec_overlaps(acb_theta_ctx_exp_zs_inv(ctx1), acb_theta_ctx_exp_zs_inv(ctx2), 6 * g))
        {
            flint_printf("FAIL (exp_zs_inv)\n");
            flint_abort();
        }
        if (!_acb_vec_overlaps(acb_theta_ctx_exp_2zs(ctx1), acb_theta_ctx_exp_2zs(ctx2), 6 * g))
        {
            flint_printf("FAIL (exp_2zs)\n");
            flint_abort();
        }
        if (!_acb_vec_overlaps(acb_theta_ctx_exp_2zs_inv(ctx1), acb_theta_ctx_exp_2zs_inv(ctx2), 6 * g))
        {
            flint_printf("FAIL (exp_2zs_inv)\n");
            _acb_vec_printd(acb_theta_ctx_exp_2zs(ctx1), 6 * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_2zs_inv(ctx1), 6 * g, 5);
            _acb_vec_printd(acb_theta_ctx_exp_2zs_inv(ctx2), 6 * g, 5);
            flint_abort();
        }
        if (!_acb_vec_overlaps(acb_theta_ctx_cs(ctx1), acb_theta_ctx_cs(ctx2), 6))
        {
            flint_printf("FAIL (cs)\n");
            flint_abort();
        }
        /* if (!_arb_vec_overlaps(acb_theta_ctx_as(ctx1), acb_theta_ctx_as(ctx2), 6 * g))
        {
            flint_printf("FAIL (as)\n");
            flint_abort();
            } */
        /* Do we need ctx_as and ctx_us ? */
        if (g > 1 && !arb_mat_overlaps(acb_theta_ctx_cho(ctx1), acb_theta_ctx_cho(ctx2)))
        {
            flint_printf("FAIL (cho)\n");
            flint_abort();
        }
        if (g > 1 && !arb_mat_overlaps(acb_theta_ctx_choinv(ctx1), acb_theta_ctx_choinv(ctx2)))
        {
            flint_printf("FAIL (choinv)\n");
            flint_abort();
        }
        if (g > 1 && !acb_mat_overlaps(acb_theta_ctx_exp_tau_inv(ctx1), acb_theta_ctx_exp_tau_inv(ctx2)))
        {
            flint_printf("FAIL (exp_tau_inv)\n");
            flint_abort();
        }
        if (g > 1 && !_arb_vec_overlaps(acb_theta_ctx_vs(ctx1), acb_theta_ctx_vs(ctx2), 6 * g))
        {
            flint_printf("FAIL (vs)\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            _arb_vec_printd(acb_theta_ctx_vs(ctx1), 6 * g, 5);
            _arb_vec_printd(acb_theta_ctx_vs(ctx2), 6 * g, 5);
            flint_abort();
        }
        if (g > 1 && !_arb_vec_overlaps(acb_theta_ctx_d0(ctx1), acb_theta_ctx_d0(ctx2), n))
        {
            flint_printf("FAIL (d0)\n");
            flint_abort();
        }
        if (g > 1 && !_arb_vec_overlaps(acb_theta_ctx_d(ctx1), acb_theta_ctx_d(ctx2), n))
        {
            flint_printf("FAIL (d)\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            _arb_vec_printd(acb_theta_ctx_d(ctx1), n, 5);
            _arb_vec_printd(acb_theta_ctx_d(ctx2), n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        acb_theta_ctx_clear(ctx1);
        acb_theta_ctx_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
