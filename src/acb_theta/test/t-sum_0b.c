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

TEST_FUNCTION_START(acb_theta_sum_0b, state)
{
    slong iter;

    /* Test: matches naive_0b */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 4);
        slong nb = 1 + n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr zs;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        acb_ptr th1, th2;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        acb_theta_ctx_tau_init(ctx_tau, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);
        th1 = _acb_vec_init(n * nb);
        th2 = _acb_vec_init(n * nb);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(zs, state, nb * g, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }
        acb_theta_sum_0b(th1, vec, nb, ctx_tau, prec);
        acb_theta_naive_0b(th2, zs, nb, tau, prec);

        if (!_acb_vec_overlaps(th1, th2, n * nb))
        {
            flint_printf("FAIL\n");
            flint_printf("\n\ng=%wd\n", g);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(zs, nb * g, 5);
            flint_printf("th1 : ");
            _acb_vec_printd(th1, n * nb, 5);
            flint_printf("th2 : ");
            _acb_vec_printd(th2, n * nb, 5);
            flint_printf("Difference: ");
            _acb_vec_sub(th1, th1, th2, n * nb, prec);
            _acb_vec_printd(th1, n * nb, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
        _acb_vec_clear(th1, n * nb);
        _acb_vec_clear(th2, n * nb);
    }

    TEST_FUNCTION_END(state);
}
