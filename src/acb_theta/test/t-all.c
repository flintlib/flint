/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_all, state)
{
    slong iter;

    /* Test: agrees with all_notransform */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n2 = 1 << (2 * g);
        slong nb = n_randint(state, 3);
        slong mprec = 50 + n_randint(state, 400);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 4);
        int sqr = n_randint(state, 2);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, test;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(nb * g);
        th = _acb_vec_init(nb * n2);
        test = _acb_vec_init(nb * n2);

        /* Sample tau not too far from reduced domain */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        acb_siegel_randtest_vec_reduced(z, state, nb, tau, 0, prec);

        /* Call theta_all at precision mprec, test against all_notransform */
        acb_theta_all(th, z, nb, tau, sqr, mprec);
        acb_theta_all_notransform(test, z, nb, tau, sqr, prec);

        if (!_acb_vec_overlaps(th, test, nb * n2))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, mprec = %wd, prec = %wd, nb = %wd, sqr = %wd, tau, z:\n",
                g, mprec, prec, nb, sqr);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, nb * g, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, nb * n2, 5);
            _acb_vec_printd(test, nb * n2, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, nb * g);
        _acb_vec_clear(th, nb * n2);
        _acb_vec_clear(test, nb * n2);
    }

    TEST_FUNCTION_END(state);
}
