/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_00, state)
{
    slong iter;

    /* Test: agrees with 00_notransform */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_randint(state, 3);
        slong prec = 100 + n_randint(state, 400);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb);
        test = _acb_vec_init(nb);

        /* Sample tau not too far from reduced domain */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        for (k = 0; k < nb; k++)
        {
            acb_siegel_randtest_vec_reduced(z + k * g, state, tau, 0, prec);
        }

        acb_theta_00(th, z, nb, tau, prec);
        acb_theta_00_notransform(test, z, nb, tau, prec);

        /*   flint_printf("g = %wd, prec = %wd, nb = %wd, tau, z:\n", g, prec, nb);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, nb * g, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nb, 5);
            _acb_vec_printd(test, nb, 5); */

        if (!_acb_vec_overlaps(th, test, nb))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, nb = %wd, tau, z:\n", g, prec, nb);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, nb * g, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nb, 5);
            _acb_vec_printd(test, nb, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, nb * g);
        _acb_vec_clear(th, nb);
        _acb_vec_clear(test, nb);
    }

    TEST_FUNCTION_END(state);
}
