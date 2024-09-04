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

TEST_FUNCTION_START(acb_theta_jet_all_notransform, state)
{
    slong iter;

    /* Test: coincides with sum_jet_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n2 = 1 << (2 * g);
        slong prec = 100 + n_randint(state, (g == 1 ? 5000 : 500));
        slong bits = n_randint(state, 4);
        slong nb = n_randint(state, 4);
        slong ord = n_randint(state, 4);
        slong nbth = acb_theta_jet_nb(ord, g);
        acb_mat_t tau;
        acb_ptr zs, th, test;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * nbth * n2);
        test = _acb_vec_init(nb * nbth * n2);
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(zs, state, nb, tau, 0, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }
        acb_theta_sum_jet_all(test, vec, nb, ctx_tau, ord, prec);

        /*flint_printf("\n\ng = %wd, prec = %wd, nb = %wd, ord = %wd\n",
            g, prec, nb, ord);
        acb_mat_printd(tau, 5);
        _acb_vec_printd(zs, nb * g, 5); */
        /*flint_printf("result of sum:\n");
          _acb_vec_printd(test, nb * n2 * nbth, 5);*/

        acb_theta_jet_all_notransform(th, zs, nb, tau, ord, prec);

        /* flint_printf("\nall_notransform: got theta:\n", res);
           _acb_vec_printd(th, nb * n2 * nbth, 5); */

        if (!_acb_vec_overlaps(th, test, nb * nbth * n2)
            || (_acb_vec_is_finite(test, nb * nbth * n2) && !_acb_vec_is_finite(th, nb * nbth * n2)))
        {
            flint_printf("FAIL\n");
            flint_printf("difference:\n");
            _acb_vec_sub(th, th, test, nb * n2 * nbth, prec);
            _acb_vec_printd(th, nb * n2 * nbth, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(th, n2 * nb * nbth);
        _acb_vec_clear(test, n2 * nb * nbth);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
    }

    TEST_FUNCTION_END(state);
}
