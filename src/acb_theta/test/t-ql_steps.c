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

TEST_FUNCTION_START(acb_theta_ql_steps, state)
{
    slong iter;

    /* Test: coincides with sum */
    for (iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 200);
        slong nb = 1 + n_randint(state, 4);
        slong nb_steps = 1 + n_randint(state, 4);
        int all = n_randint(state, 2);
        slong nbth = (all ? n * n : n);
        acb_mat_t tau;
        acb_ptr zs, t, rts, rts_all;
        acb_ptr th, th_init, test;
        arb_ptr distances;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx_z;
        slong * easy_steps;
        slong guard, hp;
        slong j;
        int res;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        t = _acb_vec_init(g);
        rts = _acb_vec_init(nb * 3 * n * nb_steps);
        if (all)
        {
            rts_all = _acb_vec_init(n * n * nb);
        }
        th = _acb_vec_init(nb * nbth);
        th_init = _acb_vec_init(3 * n * nb);
        test = _acb_vec_init(nb * nbth);
        distances = _arb_vec_init(nb * n);
        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        acb_theta_ctx_z_init(ctx_z, g);
        easy_steps = flint_malloc(nb * sizeof(slong));

        acb_siegel_randtest_compact(tau, state, 1, prec);
        acb_siegel_randtest_vec_reduced(zs + g, state, nb - 1, tau, 1, prec);
        acb_theta_eld_distances(distances, zs, nb, tau, prec);
        res = acb_theta_ql_setup(rts, rts_all, t, &guard, easy_steps, zs, nb, tau,
            distances, nb_steps, all, prec);
        hp = prec + guard;
        /* flint_printf("\n\ng = %wd, prec = %wd, nb = %wd, nb_steps = %wd, all = %wd\n", g, prec, nb, nb_steps, all);
           acb_mat_printd(tau, 5);
           _acb_vec_printd(zs, nb * g, 5); */

        if (res)
        {
            acb_theta_ql_steps_input_sum(th_init, zs, nb, t, tau, distances,
                nb_steps, easy_steps, hp);
            acb_theta_ql_steps(th, th_init, rts, rts_all, nb, nb_steps, distances,
                easy_steps, all, g, hp);

            acb_theta_ctx_tau_set(ctx_tau, tau, prec);
            for (j = 0; j < nb; j++)
            {
                acb_theta_ctx_z_set(ctx_z, zs + j * g, ctx_tau, prec);
                acb_theta_sum(test + j * nbth, ctx_z, 1, ctx_tau, distances + j * n, 1, all, 1, prec);
            }

            /* flint_printf("After %wd steps:\n", nb_steps);
               _acb_vec_printd(th, nb * nbth, 5);
               flint_printf("Result of sum:\n");
               _acb_vec_printd(test, nb * nbth, 5); */

            if (!_acb_vec_overlaps(th, test, nb * nbth))
            {
                flint_printf("FAIL\n");
                flint_printf("difference:\n");
                _acb_vec_sub(th, th, test, nb * nbth, prec);
                _acb_vec_printd(th, nb * nbth, 5);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(rts, nb * 3 * n * nb_steps);
        if (all)
        {
            _acb_vec_clear(rts_all, n * n * nb);
        }
        _acb_vec_clear(th, nbth * nb);
        _acb_vec_clear(th_init, 3 * n * nb);
        _acb_vec_clear(test, nbth * nb);
        _arb_vec_clear(distances, nb * n);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx_z);
        flint_free(easy_steps);
    }

    TEST_FUNCTION_END(state);
}
