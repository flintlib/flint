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

TEST_FUNCTION_START(acb_theta_ctx_tau_dupl, state)
{
    slong iter;

    /* Test: matches with acb_theta_ctx_tau_set with doubled input */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_theta_ctx_tau_t ctx1, ctx2;
        ulong a;

        acb_mat_init(tau, g, g);
        acb_theta_ctx_tau_init(ctx1, g);
        acb_theta_ctx_tau_init(ctx2, g);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_theta_ctx_tau_set(ctx1, tau, prec);
        acb_theta_ctx_tau_dupl(ctx1, prec);

        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        acb_theta_ctx_tau_set(ctx2, tau, prec);

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
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_4(ctx1),
                acb_theta_ctx_exp_tau_div_4(ctx2)))
        {
            flint_printf("FAIL (exp_tau_div_4)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_2(ctx1),
                acb_theta_ctx_exp_tau_div_2(ctx2)))
        {
            flint_printf("FAIL (exp_tau_div_2)\n");
            flint_abort();
        }
        if (!acb_mat_overlaps(acb_theta_ctx_exp_tau(ctx1), acb_theta_ctx_exp_tau(ctx2)))
        {
            flint_printf("FAIL (exp_tau)\n");
            flint_abort();
        }
        if (g > 1)
        {
            if (!arb_mat_overlaps(acb_theta_ctx_cho(ctx1), acb_theta_ctx_cho(ctx2)))
            {
                flint_printf("FAIL (cho)\n");
                flint_abort();
            }
            if (!arb_mat_overlaps(acb_theta_ctx_choinv(ctx1), acb_theta_ctx_choinv(ctx2)))
            {
                flint_printf("FAIL (choinv)\n");
                flint_abort();
            }
            if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_4_inv(ctx1),
                    acb_theta_ctx_exp_tau_div_4_inv(ctx2)))
            {
                flint_printf("FAIL (exp_tau_div_4_inv)\n");
                flint_abort();
            }
            if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_div_2_inv(ctx1),
                    acb_theta_ctx_exp_tau_div_2_inv(ctx2)))
            {
                flint_printf("FAIL (exp_tau_div_2_inv)\n");
                flint_abort();
            }
            if (!acb_mat_overlaps(acb_theta_ctx_exp_tau_inv(ctx1),
                    acb_theta_ctx_exp_tau_inv(ctx2)))
            {
                flint_printf("FAIL (exp_tau_inv)\n");
                flint_abort();
            }
            for (a = 0; a < n; a++)
            {
                if (!_acb_vec_overlaps(acb_theta_ctx_exp_tau_a_div_2(ctx1, a),
                        acb_theta_ctx_exp_tau_a_div_2(ctx2, a), g))
                {
                    flint_printf("FAIL (exp_tau_a_div_2)\n");
                    flint_abort();
                }
                if (!_acb_vec_overlaps(acb_theta_ctx_exp_tau_a(ctx1, a),
                        acb_theta_ctx_exp_tau_a(ctx2, a), g))
                {
                    flint_printf("FAIL (exp_tau_a)\n");
                    flint_abort();
                }
                if (!_acb_vec_overlaps(acb_theta_ctx_exp_tau_a_div_2_inv(ctx1, a),
                        acb_theta_ctx_exp_tau_a_div_2_inv(ctx2, a), g))
                {
                    flint_printf("FAIL (exp_tau_a_div_2_inv)\n");
                    flint_abort();
                }
                if (!_acb_vec_overlaps(acb_theta_ctx_exp_tau_a_inv(ctx1, a),
                        acb_theta_ctx_exp_tau_a_inv(ctx2, a), g))
                {
                    flint_printf("FAIL (exp_tau_a_inv)\n");
                    flint_abort();
                }
                if (!acb_overlaps(acb_theta_ctx_exp_a_tau_a_div_4(ctx1, a),
                        acb_theta_ctx_exp_a_tau_a_div_4(ctx2, a)))
                {
                    flint_printf("FAIL (exp_a_tau_a_div_4)\n");
                    flint_abort();
                }
            }
        }

        acb_mat_clear(tau);
        acb_theta_ctx_tau_clear(ctx1);
        acb_theta_ctx_tau_clear(ctx2);
    }

    TEST_FUNCTION_END(state);
}
