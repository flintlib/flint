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

TEST_FUNCTION_START(acb_theta_ql_all, state)
{
    slong iter;

    /* Test: agrees with naive_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        int hasz = iter % 2;
        int sqr = iter % 3;
        slong prec = (g > 1 ? 100 : 1000) + n_randint(state, 200);
        slong hprec = prec + 25;
        slong bits = 6;
        acb_mat_t tau;
        acb_ptr z, th, test;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(n * n);
        test = _acb_vec_init(n * n);

        acb_siegel_randtest_reduced(tau, state, hprec, bits);
        if (hasz)
        {
            acb_siegel_randtest_vec(z, state, g, prec);
        }

        acb_theta_ql_all(th, z, tau, sqr, prec);
        acb_theta_naive_all(test, z, 1, tau, hprec);
        if (sqr)
        {
            _acb_vec_sqr(test, test, n * n, prec);
        }

        if (!_acb_vec_overlaps(th, test, n * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, hasz = %wd, tau, z:\n",
                g, prec, hasz);
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("output:\n");
            _acb_vec_printd(th, n * n, 5);
            _acb_vec_printd(test, n * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th, n * n);
        _acb_vec_clear(test, n * n);
    }

    TEST_FUNCTION_END(state);
}

