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

TEST_FUNCTION_START(acb_theta_sum_00, state)
{
    slong iter;

    /* Test: matches sum_0b */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong mprec = 50 + n_randint(state, 200);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        slong nb = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr zs;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        acb_ptr th1, th2;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        acb_theta_ctx_tau_init(ctx_tau, 0, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);
        th1 = _acb_vec_init(nb);
        th2 = _acb_vec_init(nb * n);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(zs, state, nb, tau, 0, prec);
        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }

        /* Call sum_00 at precision mprec, test against sum_0b */
        acb_theta_sum_00(th1, vec, nb, ctx_tau, mprec);
        acb_theta_sum_0b(th2, vec, nb, ctx_tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_set(&th2[j], &th2[n * j]);
        }

        if (!_acb_vec_overlaps(th1, th2, nb))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, mprec = %wd, prec = %wd\n", g, mprec, prec);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(zs, nb * g, 5);
            flint_printf("th1 : ");
            _acb_vec_printd(th1, nb, 5);
            flint_printf("th2 : ");
            _acb_vec_printd(th2, nb, 5);
            flint_printf("Difference: ");
            _acb_vec_sub(th1, th1, th2, nb, prec);
            _acb_vec_printd(th1, nb, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
        _acb_vec_clear(th1, nb);
        _acb_vec_clear(th2, nb * n);
    }

    TEST_FUNCTION_END(state);
}
