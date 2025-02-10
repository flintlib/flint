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
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ql_jet_fd, state)
{
    slong iter;

    /* Test: coincides with sum_jet */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        int all = iter % 2;
        slong n = 1 << g;
        slong nbth = (all ? n * n : n);
        slong mprec = 100 + n_randint(state, 100);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        slong nb = n_randint(state, 3);
        slong ord = n_randint(state, 3);
        slong nbjet = acb_theta_jet_nb(ord, g);
        acb_mat_t tau;
        acb_ptr zs, th, test;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * nbth * nbjet);
        test = _acb_vec_init(nb * nbth * nbjet);
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(zs, state, nb, tau, 0, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }

        acb_theta_sum_jet(test, vec, nb, ctx_tau, ord, 1, all, prec);
        acb_theta_ql_jet_fd(th, zs, nb, tau, ord, all, mprec);

        if (!_acb_vec_overlaps(th, test, nb * nbth * nbjet)
            || !_acb_vec_is_finite(test, nb * nbth * nbjet)
            || !_acb_vec_is_finite(th, nb * nbth * nbjet))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, nb = %wd, ord = %wd, all = %wd, mprec = %wd, prec = %wd\n",
                g, nb, ord, all, mprec, prec);
            _acb_vec_printd(th, nb * nbth * nbjet, 5);
            _acb_vec_printd(test, nb * nbth * nbjet, 5);
            flint_printf("difference:\n");
            _acb_vec_sub(th, th, test, nb * nbth * nbjet, prec);
            _acb_vec_printd(th, nb * nbth * nbjet, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(th, nb * nbth * nbjet);
        _acb_vec_clear(test, nb * nbth * nbjet);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
    }

    TEST_FUNCTION_END(state);
}
