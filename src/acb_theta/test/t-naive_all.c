/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("naive_all....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with built-in genus 1 on diagonal matrices */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_pow(2, 2 * g);
        acb_mat_t tau;
        acb_mat_t tau11;
        acb_ptr z;
        slong nbz = 1 + n_randint(state, 4);
        acb_ptr th;
        acb_ptr th_test;
        acb_ptr th_g1;
        slong prec1 = 20 + n_randint(state, 200);
        slong prec = prec1; /* + n_randint(state, 200);*/
        slong mag_bits = n_randint(state, 2);
        slong k, j;
        ulong ab, a, b;

        acb_mat_init(tau, g, g);
        acb_mat_init(tau11, 1, 1);
        z = _acb_vec_init(g * nbz);
        th = _acb_vec_init(nb * nbz);
        th_test = _acb_vec_init(nb * nbz);
        th_g1 = _acb_vec_init(4 * g);

        for (k = 0; k < g; k++)
        {
            acb_siegel_randtest(tau11, state, prec, mag_bits);
            acb_set(acb_mat_entry(tau, k, k), acb_mat_entry(tau11, 0, 0));
        }
        for (k = 0; k < g * nbz; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_naive_all(th, z, nbz, tau, prec1);

        if (g == 1)
        {
            for (k = 0; k < nbz; k++)
            {
                acb_modular_theta(&th_test[4 * k + 3], &th_test[4 * k + 2],
                    &th_test[4 * k], &th_test[4 * k + 1],
                    z + k * g, acb_mat_entry(tau, 0, 0), prec);
                acb_neg(&th_test[4 * k + 3], &th_test[4 * k + 3]);
            }
        }
        else
        {
            for (j = 0; j < nbz; j++)
            {
                for (k = 0; k < g; k++)
                {
                    acb_set(acb_mat_entry(tau11, 0, 0), acb_mat_entry(tau, k, k));
                    acb_theta_naive_all(&th_g1[4 * k], &z[j * g + k], 1, tau11, prec);
                }
                /* Could use a more efficient recursive algorithm here */
                for (ab = 0; ab < n_pow(2, 2 * g); ab++)
                {
                    a = ab >> g;
                    b = ab;
                    acb_one(&th_test[j * nb + ab]);
                    for (k = g - 1; k >= 0; k--)
                    {
                        acb_mul(&th_test[j * nb + ab], &th_test[j * nb + ab],
                            &th_g1[4 * k + 2 * (a % 2) + (b % 2)], prec);
                        a = a >> 1;
                        b = b >> 1;
                    }
                }
            }
        }

        if (!_acb_vec_overlaps(th, th_test, nb * nbz))
        {
            flint_printf("FAIL: overlap\n");
            flint_printf("g = %wd, prec = %wd, nbz = %wd, tau:\n", g, prec, nbz);
            acb_mat_printd(tau, 10);
            flint_printf("z:\n");
            _acb_vec_printd(z, g * nbz, 10);
            flint_printf("th, th_test:\n");
            _acb_vec_printd(th, nb * nbz, 10);
            _acb_vec_printd(th_test, nb * nbz, 10);
            fflush(stdout);
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(th, nb * nbz);
        _acb_vec_clear(th_test, nb * nbz);
        _acb_vec_clear(th_g1, 4 * g);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
