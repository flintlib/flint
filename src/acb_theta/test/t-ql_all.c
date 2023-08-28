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

    flint_printf("ql_all....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive_all */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        int has_z = iter % 2;
        slong prec = (g > 1 ? 100 : 1000) + n_randint(state, 500);
        slong hprec = prec + 25;
        slong bits = n_randint(state, 3);
        acb_mat_t tau;
        acb_ptr z, th, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th = _acb_vec_init(n * n);
        test = _acb_vec_init(n * n);

        acb_siegel_randtest_reduced(tau, state, hprec, bits);
        if (has_z)
        {
            for (k = 0; k < g; k++)
            {
                acb_urandom(&z[k], state, hprec);
            }
        }

        acb_theta_ql_all(th, z, tau, prec);
        acb_theta_naive_all(test, z, 1, tau, hprec);

        if (!acb_is_finite(&th[0]) || !_acb_vec_overlaps(th, test, n * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, has_z = %wd, tau:\n",
                g, prec, has_z);
            acb_mat_printd(tau, 5);
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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

