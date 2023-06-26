/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("newton_sqr....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong nb = 1 << g;
        acb_mat_t tau;
        acb_ptr z;
        arf_t rad;
        slong rad_exp = -2;
        acb_ptr th2;
        acb_ptr th2_test;
        slong prec = (2 + n_randint(state, 4 - g)) * ACB_THETA_AGM_BASEPREC
            + 100 * n_randint(state, 20);
        int res;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        th2 = _acb_vec_init(2 * nb);
        th2_test = _acb_vec_init(2 * nb);
        arf_init(rad);

        acb_siegel_randtest_fund(tau, state, prec);
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, rad_exp);
        for (k = 0; k < g; k++)
        {
            acb_randtest_disk(&z[k], &z[k], rad, state, prec);
        }

        acb_theta_naive(th2_test, z, 2, tau, prec);
        acb_theta_newton_sqr(th2, z, tau, prec);

        res = 1;
        for (k = 0; k < nb; k++)
        {
            acb_sqr(&th2_test[k], &th2_test[k], prec);
            if (!acb_overlaps(&th2_test[k], &th2[k]))
                res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("th[k], th_test[k]:\n");
            for (k = 0; k < nb; k++)
            {
                acb_printd(&th2[k], 10);
                flint_printf("\n");
                acb_printd(&th2_test[k], 10);
                flint_printf("\n\n");
            }
            flint_printf("Diff:\n");
            for (k = 0; k < nb; k++)
            {
                acb_sub(&th2[k], &th2[k], &th2_test[k], prec);
                acb_printd(&th2[k], 10);
                flint_printf("\n");
            }
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        _acb_vec_clear(th2, 2 * nb);
        _acb_vec_clear(th2_test, 2 * nb);
        arf_clear(rad);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
