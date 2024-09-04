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

TEST_FUNCTION_START(acb_theta_ql_setup, state)
{
    slong iter;

    /* Test: roots are indeed nonzero */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 200);
        slong nb = 1 + n_randint(state, 4);
        slong nb_steps = n_randint(state, 6);
        int all = n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr zs, t, rts, rts_all;
        arb_ptr distances;
        slong * easy_steps;
        slong guard = 0;
        slong j, k, a;
        int res;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        t = _acb_vec_init(g);
        rts = _acb_vec_init(nb * 3 * n * nb_steps);
        rts_all = _acb_vec_init(nb * n * n);
        distances = _arb_vec_init(nb * n);
        easy_steps = flint_malloc(nb * sizeof(slong));

        acb_siegel_randtest_compact(tau, state, 1, prec);
        acb_siegel_randtest_vec_reduced(zs + g, state, nb - 1, tau, 1, prec);

        /* Compute distances */
        acb_theta_agm_distances(distances, zs, nb, tau, prec);

        res = acb_theta_ql_setup(rts, rts_all, t, &guard, easy_steps, zs, nb, tau,
            distances, nb_steps, all, prec);

        if (res)
        {
            /* Check easy_steps[0] is minimal */
            for (j = 1; j < nb; j++)
            {
                if (easy_steps[j] < easy_steps[0])
                {
                    flint_printf("FAIL (easy steps)\n");
                    flint_abort();
                }
            }

            /* Check that roots are all nonzero */
            for (j = 0; j < nb; j++)
            {
                for (k = 0; k < nb_steps; k++)
                {
                    if (all && (k == 0))
                    {
                        for (a = 0; a < n * n; a++)
                        {
                            if (acb_contains_zero(&rts_all[j * n * n + a]))
                            {
                                flint_printf("FAIL (rts_all)\n");
                                flint_printf("j = %wd, k = %wd, a = %wd, all = %wd\n", j, k, a, all);
                                acb_printd(&rts_all[j * n * n + a], 5);
                                flint_printf("\n");
                                flint_abort();
                            }
                        }
                    }
                    else
                    {
                        for (a = 0; a < n; a++)
                        {
                            if ((k < easy_steps[j])
                                && acb_contains_zero(&rts[j * (3 * n * nb_steps) + k * (3 * n) + a]))
                            {
                                flint_printf("FAIL (root 1)\n");
                                flint_abort();
                            }
                            if (k >= easy_steps[j]
                                && acb_contains_zero(&rts[j * (3 * n * nb_steps) + k * (3 * n) + n + a]))
                            {
                                flint_printf("FAIL (root 2)\n");
                                flint_printf("j = %wd, k = %wd, a = %wd\n", j, k, a);
                                acb_printd(&rts[j * (3 * n * nb_steps) + k * (3 * n) + n + a], 5);
                                flint_printf("\n");
                                flint_abort();
                            }
                            if (k >= easy_steps[j]
                                && acb_contains_zero(&rts[j * (3 * n * nb_steps) + k * (3 * n) + 2 * n + a]))
                            {
                                flint_printf("FAIL (root 3)\n");
                                flint_abort();
                            }
                        }
                    }
                }
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(rts, nb * 3 * n * nb_steps);
        _acb_vec_clear(rts_all, nb * n * n);
        _arb_vec_clear(distances, nb * n);
        flint_free(easy_steps);
    }

    TEST_FUNCTION_END(state);
}
