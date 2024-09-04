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

TEST_FUNCTION_START(acb_theta_ql_lower_dim, state)
{
    slong iter;

    /* Test: agrees with sum_a0_tilde */
    for (iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 2);
        slong n = 1 << g;
        slong s = 1 + n_randint(state, g - 1);
        slong n0 = 1 << s;
        slong nba = 1 << (g - s);
        ulong a = n_randint(state, nba);
        slong prec = 100 + n_randint(state, 100);
        int all = n_randint(state, 2);
        slong nbth = (all ? n * n : n);
        slong nbth0 = (all ? n0 * n0 : n0);
        acb_mat_t tau, tau0;
        acb_ptr z;
        arb_ptr d, d0;
        acb_ptr z0s, cofactors;
        slong * pts;
        arf_t err;
        acb_theta_ctx_tau_t ctx_tau, ctx_tau0;
        acb_theta_ctx_z_t ctx_z, ctx_z0;
        acb_ptr th, th0s, test;
        slong nb, fullprec;
        slong j;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        d = _arb_vec_init(n);
        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        acb_theta_ctx_tau_init(ctx_tau0, 1, s);
        acb_theta_ctx_z_init(ctx_z, g);
        acb_theta_ctx_z_init(ctx_z0, s);
        th = _acb_vec_init(nbth);
        test = _acb_vec_init(nbth);
        arf_init(err);

        acb_siegel_randtest_compact(tau, state, 0, prec);
        acb_siegel_randtest_vec_reduced(z, state, 1, tau, 0, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx_z, z, ctx_tau, prec);
        acb_theta_agm_distances(d, z, 1, tau, ACB_THETA_LOW_PREC);
        if (all)
        {
            acb_theta_sum_all_tilde(th, ctx_z, 1, ctx_tau, d, prec);
        }
        else
        {
            acb_theta_sum_a0_tilde(th, ctx_z, 1, ctx_tau, d, prec);
        }

        /*flint_printf("\n\ng = %wd, s = %wd, all = %wd, a = %wd, tau, z:\n", g, s, all, a);
           acb_mat_printd(tau, 5);
           _acb_vec_printd(z, g, 5); */

        res = acb_theta_ql_lower_dim(&z0s, &cofactors, &pts, &nb, err,
            &fullprec, z, tau, d, s, a, prec);

        /*flint_printf("ql_lower_dim: got res = %wd, nb = %wd, fullprec = %wd\n", res, nb, fullprec); */

        th0s = _acb_vec_init(nb * nbth0);
        d0 = _arb_vec_init(nb * n0);
        acb_mat_window_init(tau0, tau, 0, 0, s, s);

        acb_theta_ctx_tau_set(ctx_tau0, tau0, prec);
        acb_theta_agm_distances(d0, z0s, nb, tau0, ACB_THETA_LOW_PREC);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(ctx_z0, z0s + j * s, ctx_tau0, prec);
            if (all)
            {
                acb_theta_sum_all_tilde(th0s + j * nbth0, ctx_z0, 1, ctx_tau0, d0 + j * n0, prec);
            }
            else
            {
                acb_theta_sum_a0_tilde(th0s + j * nbth0, ctx_z0, 1, ctx_tau0, d0 + j * n0, prec);
            }
        }

        if (res)
        {
            acb_theta_ql_recombine(test, th0s, cofactors, pts, nb, err, fullprec,
                s, a, all, g, prec);
            /*_acb_vec_printd(th, nbth, 5);
              _acb_vec_printd(test, nbth, 5); */
            if (all)
            {
                for (j = 0; j < n * n; j++)
                {
                    if ((j >> g) % (1 << (g - s)) == a
                        && !acb_overlaps(&test[j], &th[j]))
                    {
                        flint_printf("FAIL (all, ab = %wd)\n", j);
                        flint_abort();
                    }
                }
            }
            else
            {
                for (j = 0; j < n; j++)
                {
                    if (j % (1 << (g - s)) == a
                        && !acb_overlaps(&test[j], &th[j]))
                    {
                        flint_printf("FAIL (a = %wd)\n", j);
                        flint_abort();
                    }
                }
            }
        }

        acb_mat_clear(tau);
        acb_mat_window_clear(tau0);
        _acb_vec_clear(z, g);
        _arb_vec_clear(d, n);
        _arb_vec_clear(d0, nb * n0);
        _acb_vec_clear(z0s, nb * s);
        _acb_vec_clear(cofactors, nb);
        arf_clear(err);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_tau_clear(ctx_tau0);
        acb_theta_ctx_z_clear(ctx_z);
        acb_theta_ctx_z_clear(ctx_z0);
        _acb_vec_clear(th, nbth);
        _acb_vec_clear(test, nbth);
        _acb_vec_clear(th0s, nb * nbth0);
        flint_free(pts);
    }

    TEST_FUNCTION_END(state);
}
