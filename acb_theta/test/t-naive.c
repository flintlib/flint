/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"
#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with genus 1; duplication formula */
    for (iter = 0; iter < 50 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_pow(2, g);
        acb_mat_t tau;
        acb_ptr z;
        arf_t rad;
        slong rad_exp = 0;
        acb_ptr th;
        acb_ptr th_dupl;
        acb_ptr th_test;
        slong prec = 20 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        int res;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        arf_init(rad);
        th = _acb_vec_init(2 * nb);
        th_dupl = _acb_vec_init(2 * nb * nb);
        th_test = _acb_vec_init(2 * nb * nb);

        acb_siegel_randtest(tau, state, prec, mag_bits);
        arf_one(rad);
        arf_mul_2exp_si(rad, rad, rad_exp);
        for (k = 0; k < g; k++)
        {
            acb_randtest_disk(&z[k], &z[k], rad, state, prec);
        }

        acb_theta_naive(th, z, 2, tau, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        acb_theta_naive_all(th_test, z, 2, tau, prec);

        if (g == 1)
        {
            acb_modular_theta(&th_dupl[3], &th_dupl[2],
                              &th_dupl[0], &th_dupl[1], z, acb_mat_entry(tau,
                                                                         0, 0),
                              prec);
            acb_neg(&th_dupl[3], &th_dupl[3]);
            acb_theta_naive(th, z, 1, tau, prec);
        }
        else
        {
            acb_theta_dupl_all(th_dupl, th, g, prec);
            for (k = 0; k < 2 * nb * nb; k++)
            {
                acb_sqr(&th_test[k], &th_test[k], prec);
            }
            acb_theta_naive(th, z, 1, tau, prec);
            for (k = 0; k < nb; k++)
            {
                acb_sqr(&th[k], &th[k], prec);
            }
        }

        res = 1;
        for (k = 0; k < nb * nb; k++)
        {
            if (!acb_overlaps(&th_dupl[k], &th_test[k]))
            {
                flint_printf("no overlap at k=%wd\n", k);
                res = 0;
            }
        }
        if (!res)
        {
            flint_printf("FAIL: overlap\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            for (k = 0; k < g; k++)
            {
                acb_printd(&z[k], 10);
                flint_printf("\n");
            }
            flint_printf("th_dupl[k], th_test[k], th[k]:\n");
            for (k = 0; k < 2 * nb * nb; k++)
            {
                acb_printd(&th_dupl[k], 50);
                flint_printf("\n");
                acb_printd(&th_test[k], 50);
                flint_printf("\n");
                if (k < nb)
                {
                    acb_printd(&th[k], 10);
                    flint_printf("\n");
                }
                flint_printf("\n");
            }
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        arf_clear(rad);
        _acb_vec_clear(th, 2 * nb);
        _acb_vec_clear(th_dupl, 2 * nb * nb);
        _acb_vec_clear(th_test, 2 * nb * nb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
