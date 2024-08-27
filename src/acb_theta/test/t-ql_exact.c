/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ql_exact, state)
{
    slong iter;

    /* Test: coincides with sum_a0_tilde */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 500);
        slong bits = n_randint(state, 4);
        slong nb = 1 + n_randint(state, 4);
        slong * pattern;
        int all = n_randint(state, 2);
        int shifted_prec = n_randint(state, 2);
        slong nbth = (all ? n * n : n);
        acb_mat_t tau;
        arb_t y;
        acb_ptr zs, th, test;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        arb_ptr distances;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        arb_init(y);
        zs = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * nbth);
        test = _acb_vec_init(nb * nbth);
        acb_theta_ctx_tau_init(ctx_tau, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);
        distances = _arb_vec_init(nb * n);
        pattern = flint_malloc(g * sizeof(slong));

        /* Sample tau with reasonable imaginary part */
        res = 0;
        while(!res)
        {
            acb_siegel_randtest_reduced(tau, state, prec, bits);
            arb_sub_si(y, acb_imagref(acb_mat_entry(tau, g - 1, g - 1)), 200, prec);
            res = arb_is_negative(y);
        }
        acb_siegel_randtest_vec(zs + g, state, (nb - 1) * g, prec);
        for (j = 0; j < g; j++)
        {
            pattern[j] = n_randint(state, 4);
        }

        /* Strip tau, z of errors */
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                acb_get_mid(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k));
                acb_set(acb_mat_entry(tau, k, j), acb_mat_entry(tau, j, k));
            }
        }
        for (j = 0; j < nb * g; j++)
        {
            acb_get_mid(&zs[j], &zs[j]);
        }

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
            if (shifted_prec)
            {
                acb_theta_dist_a0(distances + j * n, zs + j * g, tau, prec);
                if (all)
                {
                    acb_theta_sum_all_tilde(test + j * n * n, &vec[j], 1, ctx_tau, distances + j * n, prec);
                }
                else
                {
                    acb_theta_sum_a0_tilde(test + j * n, &vec[j], 1, ctx_tau, distances + j * n, prec);
                }
            }
        }
        if (!shifted_prec)
        {
            /* distances are set to zero */
            if (all)
            {
                acb_theta_sum_all_tilde(test, vec, nb, ctx_tau, distances, prec);
            }
            else
            {
                acb_theta_sum_a0_tilde(test, vec, nb, ctx_tau, distances, prec);
            }
        }

        /* flint_printf("\n\ng = %wd, prec = %wd, nb = %wd, all = %wd, shifted_prec = %wd\n",
            g, prec, nb, all, shifted_prec);
        acb_mat_printd(tau, 5);
        _acb_vec_printd(zs, nb * g, 5);
        flint_printf("result of sum:\n");
        _acb_vec_printd(test, nbth * nb, 5);
        flint_printf("pattern:\n");
        for (j = 0; j < g; j++)
        {
            flint_printf("%wd -> %wd\n", j, pattern[j]);
            } */

        res = acb_theta_ql_exact(th, zs, nb, tau, pattern, all, shifted_prec, prec);

        /* flint_printf("\nresult of ql_exact: %wd, got theta:\n", res);
           _acb_vec_printd(th, nbth * nb, 5); */

        if (!res
            || !_acb_vec_overlaps(th, test, nb * nbth)
            || (_acb_vec_is_finite(test, nb * nbth) && !_acb_vec_is_finite(th, nb * nbth)))
        {
            flint_printf("FAIL\n");
            flint_printf("difference:\n");
            _acb_vec_sub(th, th, test, nb * nbth, prec);
            _acb_vec_printd(th, nb * nbth, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        arb_clear(y);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(th, nbth * nb);
        _acb_vec_clear(test, nbth * nb);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
        _arb_vec_clear(distances, nb * n);
        flint_free(pattern);
    }

    TEST_FUNCTION_END(state);
}
