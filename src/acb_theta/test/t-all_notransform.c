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

TEST_FUNCTION_START(acb_theta_all_notransform, state)
{
    slong iter;

    /* Test: coincides with sum_all_tilde */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 500);
        slong nb = 1 + n_randint(state, 4);
        int sqr = n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr zs, th, test;
        arb_ptr distances;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        slong j;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * n * n);
        test = _acb_vec_init(nb * n * n);
        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);
        distances = _arb_vec_init(n);

        /* Sample tau with reasonable imaginary part */
        acb_siegel_randtest_compact(tau, state, 0, prec);
        acb_siegel_randtest_vec_reduced(zs + g, state, nb - 1, tau, 0, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
        }
        acb_theta_sum_all_tilde(test, vec, nb, ctx_tau, distances, prec);
        for (j = 0; j < nb; j++)
        {
            _acb_vec_scalar_mul_arb(test + j * n * n, test + j * n * n, n * n,
                &(&vec[j])->u, prec);
        }
        if (sqr)
        {
            for (j = 0; j < n * n * nb; j++)
            {
                acb_sqr(&test[j], &test[j], prec);
            }
        }

        /*flint_printf("\n\ng = %wd, prec = %wd, nb = %wd, sqr = %wd\n",
            g, prec, nb, sqr);
         acb_mat_printd(tau, 5);
        _acb_vec_printd(zs, nb * g, 5);
        flint_printf("result of sum:\n");
        _acb_vec_printd(test, n * n * nb, 5); */

        acb_theta_all_notransform(th, zs, nb, tau, sqr, prec);

        /* flint_printf("\nall_notransform: got theta:\n", res);
           _acb_vec_printd(th, n * n * nb, 5); */

        if (!_acb_vec_overlaps(th, test, nb * n * n)
            || (_acb_vec_is_finite(test, nb * n * n) && !_acb_vec_is_finite(th, nb * n * n)))
        {
            flint_printf("FAIL\n");
            flint_printf("difference:\n");
            _acb_vec_sub(th, th, test, nb * n * n, prec);
            _acb_vec_printd(th, nb * n * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(th, n * n * nb);
        _acb_vec_clear(test, n * n * nb);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
        _arb_vec_clear(distances, n);
    }

    TEST_FUNCTION_END(state);
}
