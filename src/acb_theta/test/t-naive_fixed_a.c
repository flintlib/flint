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

TEST_FUNCTION_START(acb_theta_naive_fixed_a, state)
{
    slong iter;

    /* Test: agrees with naive_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong nbz = 1 + n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, th_all, th_test;
        slong prec = 20 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 2);
        slong k, a;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g * nbz);
        th = _acb_vec_init(n * nbz);
        th_all = _acb_vec_init(n * n * nbz);
        th_test = _acb_vec_init(n * nbz);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g * nbz, prec);

        acb_theta_naive_all(th_all, z, nbz, tau, prec);

        for (a = 0; a < n; a++)
        {
            acb_theta_naive_fixed_a(th, a, z, nbz, tau, prec);
            for (k = 0; k < nbz; k++)
            {
                _acb_vec_set(th_test + k * n, th_all + k * n * n + a * n, n);
            }
            if (!_acb_vec_overlaps(th, th_test, n * nbz))
            {
                flint_printf("FAIL\n");
                flint_printf("g = %wd, prec = %wd, nbz = %wd, a = %wd, tau:\n",
                    g, prec, nbz, a);
                acb_mat_printd(tau, 5);
                flint_printf("z:\n");
                _acb_vec_printd(z, g * nbz, 10);
                flint_printf("th, th_test:\n");
                _acb_vec_printd(th, n * nbz, 10);
                _acb_vec_printd(th_test, n * nbz, 10);
                flint_printf("th_all:\n");
                _acb_vec_printd(th_all, n * n * nbz, 10);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(th, n * nbz);
        _acb_vec_clear(th_all, n * n * nbz);
        _acb_vec_clear(th_test, n * nbz);
    }

    TEST_FUNCTION_END(state);
}
