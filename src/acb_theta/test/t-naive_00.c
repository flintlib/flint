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

TEST_FUNCTION_START(acb_theta_naive_00, state)
{
    slong iter;

    /* Test: agrees with first entry of naive_0b */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        acb_mat_t tau;
        acb_ptr z;
        slong nbz = 1 + n_randint(state, 4);
        acb_ptr th, th_0b, test;
        slong prec1 = 100 + n_randint(state, 200);
        slong prec = prec1 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 10);
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g * nbz);
        th = _acb_vec_init(nbz);
        th_0b = _acb_vec_init(n * nbz);
        test = _acb_vec_init(nbz);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g * nbz, prec);

        acb_theta_naive_00(th, z, nbz, tau, prec1);
        acb_theta_naive_0b(th_0b, z, nbz, tau, prec);
        for (k = 0; k < nbz; k++)
        {
            acb_set(&test[k], &th_0b[k * n]);
        }

        if (!_acb_vec_overlaps(th, test, nbz))
        {
            flint_printf("FAIL: overlap\n");
            flint_printf("g = %wd, prec1 = %wd, prec = %wd, nbz = %wd, tau:\n",
                g, prec1, prec, nbz);
            acb_mat_printd(tau, 5);
            flint_printf("z:\n");
            _acb_vec_printd(z, g * nbz, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nbz, 5);
            _acb_vec_printd(test, nbz, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(th, nbz);
        _acb_vec_clear(th_0b, n * nbz);
        _acb_vec_clear(test, nbz);
    }

    TEST_FUNCTION_END(state);
}
