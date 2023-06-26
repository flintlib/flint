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

    flint_printf("dupl_z....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 5 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << (2 * g);
        slong prec = 200 + n_randint(state, 1000);

        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th;
        acb_ptr dupl;
        acb_ptr test;
        arf_t rad;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        th = _acb_vec_init(2 * n);
        dupl = _acb_vec_init(2 * n);
        test = _acb_vec_init(2 * n);
        arf_init(rad);

        acb_siegel_randtest_fund(tau, state, prec);
        arf_one(rad);
        for (k = 0; k < g; k++)
        {
            acb_randtest_disk(&z[k], &z[k], rad, state, prec);
        }

        acb_theta_naive_all(th, z, 2, tau, prec);
        acb_theta_dupl_z(dupl, th, g, prec);

        _acb_vec_scalar_mul_2exp_si(z, z, g, 1);
        acb_theta_naive_all(test, z, 2, tau, prec);

        res = 1;
        for (k = 0; k < 2 * n; k++)
        {
            if (!acb_overlaps(&dupl[k], &test[k]))
                res = 0;
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau:\n");
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            for (k = 0; k < g; k++)
            {
                acb_printd(&z[k], 10);
                flint_printf("\n");
            }
            flint_printf("Before dupl:\n");
            for (k = 0; k < 2 * n; k++)
            {
                acb_printd(&th[k], 10);
                flint_printf("\n");
            }
            flint_printf("Comparison:\n");
            for (k = 0; k < 2 * n; k++)
            {
                acb_printd(&test[k], 10);
                flint_printf("\n");
                acb_printd(&dupl[k], 10);
                flint_printf("\n");
            }

            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        _acb_vec_clear(th, 2 * n);
        _acb_vec_clear(dupl, 2 * n);
        _acb_vec_clear(test, 2 * n);
        arf_clear(rad);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
