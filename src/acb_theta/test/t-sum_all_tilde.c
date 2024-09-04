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
#include "acb_modular.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_sum_all_tilde, state)
{
    slong iter;

    /* Test: matches acb_modular_theta on diagonal matrices */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong n2 = 1 << (2 * g);
        slong mprec = 100 + n_randint(state, 100);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        acb_mat_t tau, tau11;
        acb_ptr z;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx;
        arb_ptr d;
        acb_ptr th, test, aux;
        slong j, ab, a1b1;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau11, 1, 1);
        z = _acb_vec_init(g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        acb_theta_ctx_z_init(ctx, g);
        d = _arb_vec_init(n);
        th = _acb_vec_init(n2);
        test = _acb_vec_init(n2);
        aux = _acb_vec_init(4);

        for (j = 0; j < g; j++)
        {
            acb_siegel_randtest_reduced(tau11, state, prec, bits);
            acb_set(acb_mat_entry(tau, j, j), acb_mat_entry(tau11, 0, 0));
        }
        acb_theta_ctx_tau_set(ctx_tau, tau, prec);

        acb_siegel_randtest_vec_reduced(z, state, tau, 0, prec);
        acb_theta_ctx_z_set(ctx, z, ctx_tau, prec);
        acb_theta_dist_a0(d, z, tau, prec);

        /* Call sum_all_tilde at precision mprec, test against modular_theta */
        acb_theta_sum_all_tilde(th, ctx, 1, ctx_tau, d, mprec);
        _acb_vec_scalar_mul_arb(th, th, n2, acb_theta_ctx_u(ctx), prec);

        for (j = 0; j < n2; j++)
        {
            acb_one(&test[j]);
        }
        for (j = 0; j < g; j++)
        {
            acb_modular_theta(&aux[3], &aux[2], &aux[0], &aux[1],
                &z[j], acb_mat_entry(tau, j, j), prec);
            acb_neg(&aux[3], &aux[3]);
            for (ab = 0; ab < n2; ab++)
            {
                a1b1 = 2 * ((ab >> (2 * g - j - 1)) % 2) + ((ab >> (g - j - 1)) % 2);
                acb_mul(&test[ab], &test[ab], &aux[a1b1], prec);
            }
        }

        if (!_acb_vec_overlaps(th, test, n2))
        {
            flint_printf("FAIL\n");
            flint_printf("g=%wd\n", g);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("th: ");
            _acb_vec_printd(th, n2, 5);
            flint_printf("test: ");
            _acb_vec_printd(test, n2, 5);
            flint_printf("Difference: ");
            _acb_vec_sub(th, th, test, n2, prec);
            _acb_vec_printd(th, n2, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(z, g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx);
        _arb_vec_clear(d, n);
        _acb_vec_clear(th, n2);
        _acb_vec_clear(test, n2);
        _acb_vec_clear(aux, 4);
    }

    TEST_FUNCTION_END(state);
}
