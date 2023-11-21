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

TEST_FUNCTION_START(acb_theta_naive_fixed_ab, state)
{
    slong iter;

    /* Test: agrees with naive_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_pow(2, g);
        slong nbz = 1 + n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, th_all, th_test;
        slong prec = 20 + n_randint(state, 100);
        slong mag_bits = n_randint(state, 2);
        ulong ab;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g * nbz);
        th = _acb_vec_init(nbz);
        th_all = _acb_vec_init(nb * nb * nbz);
        th_test = _acb_vec_init(nbz);

        acb_siegel_randtest_reduced(tau, state, prec, mag_bits);
        acb_siegel_randtest_vec(z, state, g * nbz, prec);

        acb_theta_naive_all(th_all, z, nbz, tau, prec);

        for (ab = 0; ab < nb * nb; ab++)
        {
            acb_theta_naive_fixed_ab(th, ab, z, nbz, tau, prec);
            for (k = 0; k < nbz; k++)
            {
                acb_set(&th_test[k], &th_all[k * nb * nb + ab]);
            }
            if (!_acb_vec_overlaps(th, th_test, nbz))
            {
                flint_printf("FAIL\n");
                flint_printf("g = %wd, prec = %wd, nbz = %wd, ab = %wd, tau:\n",
                    g, prec, nbz, ab);
                acb_mat_printd(tau, 5);
                flint_printf("z:\n");
                _acb_vec_printd(z, g * nbz, 10);
                flint_printf("th, th_test:\n");
                _acb_vec_printd(th, nbz, 10);
                _acb_vec_printd(th_test, nbz, 10);
                flint_printf("th_all:\n");
                _acb_vec_printd(th_all, nb * nb * nbz, 10);
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(th, nbz);
        _acb_vec_clear(th_all, nb * nb * nbz);
        _acb_vec_clear(th_test, nbz);
    }

    TEST_FUNCTION_END(state);
}
