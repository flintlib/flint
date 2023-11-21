/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_naive_all, state)
{
    slong iter;

    /* Test: agrees with built-in genus 1 on diagonal matrices */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 2);
        slong nb = n_pow(2, 2 * g);
        slong prec1 = ACB_THETA_LOW_PREC + n_randint(state, 200);
        slong prec = prec1 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 2);
        slong nbz = 1 + n_randint(state, 4);
        acb_mat_t tau, tau11;
        acb_ptr z, th, th_test, th_g1;
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
            acb_siegel_randtest_reduced(tau11, state, prec, mag_bits);
            acb_set(acb_mat_entry(tau, k, k), acb_mat_entry(tau11, 0, 0));
        }
        acb_siegel_randtest_vec(z, state, g * nbz, prec);

        acb_theta_naive_all(th, z, nbz, tau, prec1);

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
            flint_printf("difference:\n");
            _acb_vec_sub(th_test, th_test, th, nb * nbz, prec);
            _acb_vec_printd(th_test, nb * nbz, 10);
            flint_abort();
        }

        acb_mat_clear(tau);
        acb_mat_clear(tau11);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(th, nb * nbz);
        _acb_vec_clear(th_test, nb * nbz);
        _acb_vec_clear(th_g1, 4 * g);
    }

    TEST_FUNCTION_END(state);
}
