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

TEST_FUNCTION_START(acb_theta_dist_a0, state)
{
    slong iter;

    /* Test: find zero value when z = tau a/2 + real stuff */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong prec = ACB_THETA_LOW_PREC;
        slong hprec = 200;
        slong bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z;
        arb_ptr d;
        arb_t c;
        ulong a = n_randint(state, n);
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        d = _arb_vec_init(n);
        arb_init(c);

        acb_siegel_randtest_reduced(tau, state, hprec, bits);
        acb_theta_char_get_acb(z, a, g);
        acb_mat_vector_mul_col(z, tau, z, prec);
        for (k = 0; k < g; k++)
        {
            arb_urandom(c, state, prec);
            acb_add_arb(&z[k], &z[k], c, prec);
        }

        acb_theta_dist_a0(d, z, tau, prec);

        if (!arb_contains_zero(&d[a]))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, a = %wd, tau:\n", g, a);
            acb_mat_printd(tau, 5);
            flint_printf("distances:\n");
            _arb_vec_printd(d, n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _arb_vec_clear(d, n);
        arb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
