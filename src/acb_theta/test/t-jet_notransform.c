/*
    Copyright (C) 2025 Jean Kieffer

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

TEST_FUNCTION_START(acb_theta_jet_notransform, state)
{
    slong iter;

    /* Test: coincides with sum_jet */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        int all = iter % 2;
        ulong ab = n_randint(state, n);
        slong nbth = (all ? n * n : 1);
        slong mprec = 50 + n_randint(state, 100);
        slong prec = mprec + 25;
        slong bits = n_randint(state, 4);
        slong nb = n_randint(state, 3);
        slong ord = n_randint(state, 3);
        slong nbjet = acb_theta_jet_nb(ord, g);
        int sqr = n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr zs, th, test;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * nbth * nbjet);
        test = _acb_vec_init(nb * n * n * nbjet);
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(zs, state, nb, tau, 0, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }

        acb_theta_jet_notransform(th, zs, nb, tau, ord, ab, all, sqr, mprec);
        acb_theta_sum_jet(test, vec, nb, ctx_tau, ord, 1, 1, prec);
        if (!all)
        {
            for (j = 0; j < nb; j++)
            {
                _acb_vec_set(test + j * nbjet, test + j * n * n * nbjet + ab * nbjet, nbjet);
            }
        }
        if (ord == 0 && sqr)
        {
            _acb_vec_sqr(test, test, nb * nbth, prec);
        }

        if (!_acb_vec_overlaps(th, test, nb * nbth * nbjet)
            || (_acb_vec_is_finite(test, nb * nbth * nbjet)
                && !_acb_vec_is_finite(th, nb * nbth * nbjet)))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, nb = %wd, ord = %wd, all = %wd, ab = %wd, sqr = %wd, mprec = %wd, prec = %wd\n",
            g, nb, ord, all, ab, sqr, mprec, prec);
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
        _acb_vec_clear(test, nb * n * n * nbjet);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
    }

    TEST_FUNCTION_END(state);
}
