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

TEST_FUNCTION_START(acb_theta_sum_a0, state)
{
    slong iter;

    /* Test: matches naive_fixed_ab */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 4);
        int z_is_real = iter % 2;
        acb_theta_ctx_t ctx;
        acb_mat_t tau;
        acb_ptr z, t, all_zs;
        acb_ptr th1, th2, res;
        slong j, k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        all_zs = _acb_vec_init(6 * g);
        acb_theta_ctx_init(ctx, 6, g);
        th1 = _acb_vec_init(6 * n);
        th2 = _acb_vec_init(6 * n);
        res = _acb_vec_init(6);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g, prec);
        if (z_is_real)
        {
            for (j = 0; j < g; j++)
            {
                arb_zero(acb_imagref(&z[j]));
            }
        }
        for (j = 0; j < g; j++)
        {
            arb_urandom(acb_realref(&t[j]), state, prec);
        }

        acb_theta_ctx_set_tau(ctx, tau, prec);
        acb_theta_ctx_set_z_ql(ctx, z, prec);
        acb_theta_ctx_set_t(ctx, t, prec);

        if (z_is_real)
        {
            acb_theta_sum_a0(th1, ctx, 0, 6, 1, prec);
        }
        else
        {
            acb_theta_sum_a0(th1, ctx, 0, 3, 1, prec);
            acb_theta_sum_a0(th1 + 3 * n, ctx, 3, 3, 0, prec);
        }

        _acb_vec_set(all_zs + g, t, g);
        _acb_vec_scalar_mul_2exp_si(all_zs + 2 * g, t, g, 1);
        for (j = 0; j < 3; j++)
        {
            _acb_vec_add(all_zs + (3 + j) * g, all_zs + j * g, z, g, prec);
        }
        for (j = 0; j < n; j++)
        {
            acb_theta_naive_fixed_ab(res, j << g, all_zs, 6, tau, prec);
            for (k = 0; k < 6; k++)
            {
                acb_set(&th2[k * n + j], &res[k]);
            }
        }

        if (!_acb_vec_overlaps(th1, th2, 6 * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g=%wd\n", g);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(all_zs, 6 * g, 5);
            flint_printf("th1: ");
            _acb_vec_printd(th1, 6 * n, 5);
            flint_printf("th2: ");
            _acb_vec_printd(th2, 6 * n, 5);
            flint_printf("Difference: ");
            _acb_vec_sub(th1, th1, th2, 6 * n, prec);
            _acb_vec_printd(th1, 6 * n, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(all_zs, 6 * g);
        acb_theta_ctx_clear(ctx);
        _acb_vec_clear(th1, 6 * n);
        _acb_vec_clear(th2, 6 * n);
        _acb_vec_clear(res, 6);
    }

    TEST_FUNCTION_END(state);
}
