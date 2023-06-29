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

    flint_printf("naive_all....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with built-in genus 1 on diagonal matrices */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_pow(2, 2 * g);
        acb_mat_t tau;
        acb_mat_t tau11;
        acb_ptr z;
        acb_ptr th;
        acb_ptr th_test;
        acb_ptr th_g1;
        slong prec = 20 + n_randint(state, 500);
        slong mag_bits = n_randint(state, 2);
        slong k;
        ulong ab, a, b;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau11, 1, 1);
        z = _acb_vec_init(g);
        th = _acb_vec_init(nb);
        th_test = _acb_vec_init(nb);
        th_g1 = _acb_vec_init(4 * g);

        for (k = 0; k < g; k++)
        {
            acb_siegel_randtest(tau11, state, prec, mag_bits);
            acb_set(acb_mat_entry(tau, k, k), acb_mat_entry(tau11, 0, 0));
        }
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_naive_all(th, z, 1, tau, prec);

        if (g == 1)
        {
            acb_modular_theta(&th_test[3], &th_test[2], &th_test[0], &th_test[1],
                z, acb_mat_entry(tau, 0, 0), prec);
            acb_neg(&th_test[3], &th_test[3]);
        }
        else
        {
            for (k = 0; k < g; k++)
            {
                acb_set(acb_mat_entry(tau11, 0, 0), acb_mat_entry(tau, k, k));
                acb_theta_naive_all(&th_g1[4 * k], &z[k], 1, tau11, prec);
            }
            /* Could use a more efficient recursive algorithm here */
            for (ab = 0; ab < n_pow(2, 2 * g); ab++)
            {
                a = ab >> g;
                b = ab;
                acb_one(&th_test[ab]);
                for (k = g - 1; k >= 0; k--)
                {
                    acb_mul(&th_test[ab], &th_test[ab],
                        &th_g1[4 * k + 2 * (a % 2) + (b % 2)], prec);
                    a = a >> 1;
                    b = b >> 1;
                }
            }
        }
        
        if (!_acb_vec_overlaps(th, th_test, nb))
        {
            flint_printf("FAIL: overlap\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            _acb_vec_printd(z, g, 10);
            flint_printf("th, th_test:\n");
            _acb_vec_printd(th, nb, 10);
            _acb_vec_printd(th_test, nb, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, nb);
        _acb_vec_clear(th_test, nb);
        _acb_vec_clear(th_g1, 4 * g);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
