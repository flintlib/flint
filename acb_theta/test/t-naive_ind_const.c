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

    flint_printf("naive_ind_const....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with genus 1 theta */
    for (iter = 0; iter < 200 * arb_test_multiplier(); iter++)
    {
        slong g = 1;
        acb_mat_t tau;
        acb_t t;
        arb_t eps;
        acb_ptr th;
        acb_ptr th_test;
        acb_t z;
        ulong ab;
        slong prec = 20 + n_randint(state, 2000);
        slong mag_bits = n_randint(state, 10);
        int res;
        slong k;

        acb_mat_init(tau, g, g);
        th = _acb_vec_init(4);
        th_test = _acb_vec_init(4);
        acb_init(t);
        acb_init(z);
        arb_init(eps);

        arb_one(eps);
        arb_mul_2exp_si(eps, eps, -n_randint(state, 10));   /* Not too small */

        arb_randtest_precise(acb_realref(t), state, prec, mag_bits);
        arb_randtest_precise(acb_imagref(t), state, prec, mag_bits);
        arb_sqr(acb_imagref(t), acb_imagref(t), prec);
        arb_add(acb_imagref(t), acb_imagref(t), eps, prec);
        acb_set(acb_mat_entry(tau, 0, 0), t);

        acb_zero(z);
        acb_modular_theta(&th_test[3], &th_test[2],
                          &th_test[0], &th_test[1], z, t, prec);

        for (ab = 0; ab < 4; ab++)
        {
            acb_theta_naive_ind_const(&th[ab], ab, tau, prec);
        }
        res = 1;
        for (k = 0; k < 4; k++)
        {
            if (!acb_overlaps(&th[k], &th_test[k]))
                res = 0;
        }

        if (!res)
        {
            flint_printf("FAIL: no overlap\n");
            flint_printf("prec = %wd, tau:\n", prec);
            acb_mat_printd(tau, 10);
            flint_printf("\n");
            flint_printf("th_test[k], th[k] for k = 0 to 3:\n");
            for (k = 0; k < 4; k++)
            {
                acb_printd(&th_test[k], 30);
                flint_printf("\n");
                acb_printd(&th[k], 30);
                flint_printf("\n");
            }
        }

        acb_mat_clear(tau);
        acb_clear(t);
        _acb_vec_clear(th, 4);
        _acb_vec_clear(th_test, 4);
        acb_clear(z);
        arb_clear(eps);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
