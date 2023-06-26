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

    flint_printf("all_sqr....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithm */
    for (iter = 0; iter < 5 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << (2 * g);
        slong prec = 2000 + n_randint(state, 2000);
        slong mag_bits = 1 + n_randint(state, 3);

        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th2;
        acb_ptr test;
        arf_t rad;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        th2 = _acb_vec_init(2 * n);
        test = _acb_vec_init(2 * n);
        arf_init(rad);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, 2);
        for (k = 0; k < g; k++)
        {
            acb_randtest_disk(&z[k], &z[k], rad, state, prec);
        }

        /* Force unbalancedness */
        for (j = 0; j < g; j++)
        {
            for (k = 0; k < g; k++)
            {
                acb_mul_2exp_si(acb_mat_entry(tau, j, k),
                                acb_mat_entry(tau, j, k), j + k);
            }
        }

        flint_printf("\nNew matrix:\n");
        acb_mat_printd(tau, 10);
        flint_printf("New z:\n");
        for (k = 0; k < g; k++)
        {
            acb_printd(&z[k], 10);
            flint_printf("\n");
        }

        acb_theta_all_sqr(th2, z, tau, prec);
        acb_theta_naive_all(test, z, 2, tau, prec);
        acb_theta_vecsqr(test, test, 2 * n, prec);

        flint_printf("Final result:\n");
        for (k = 0; k < 2 * n; k++)
        {
            acb_printd(&th2[k], 10);
            flint_printf("\n");
            acb_printd(&test[k], 10);
            flint_printf("\n");
        }

        res = 1;
        for (k = 0; k < 2 * n; k++)
        {
            if (!acb_overlaps(&th2[k], &test[k]))
            {
                flint_printf("No overlap at k=%wd, difference\n", k);
                acb_sub(&th2[k], &th2[k], &test[k], prec);
                acb_printd(&th2[k], 10);
                flint_printf("\n");
                res = 0;
            }
        }
        if (!res)
        {
            flint_printf("FAIL (overlap)\n");
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        _acb_vec_clear(th2, 2 * n);
        _acb_vec_clear(test, 2 * n);
        arf_clear(rad);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
