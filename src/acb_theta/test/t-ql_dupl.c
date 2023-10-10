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

    flint_printf("ql_dupl....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive_all */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong lp = ACB_THETA_LOW_PREC;
        slong prec = 100;
        acb_mat_t tau;
        acb_ptr z, th0, thz, th2, test;
        arb_ptr d0, d;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th0 = _acb_vec_init(n);
        thz = _acb_vec_init(n);
        th2 = _acb_vec_init(n * n);
        test = _acb_vec_init(n * n);
        d0 = _arb_vec_init(n);
        d = _arb_vec_init(n);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);

        acb_theta_dist_a0(d0, z, tau, lp);
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_fixed_ab(&th0[k], k << g, z, 1, tau, prec);
        }

        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_dist_a0(d, z, tau, lp);
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_fixed_ab(&thz[k], k << g, z, 1, tau, prec);
        }

        acb_theta_ql_dupl(th2, th0, thz, d0, d, g, prec);

        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        _acb_vec_scalar_mul_2exp_si(z, z, g, -1);
        acb_theta_naive_all(test, z, 1, tau, prec);
        _acb_vec_sqr(test, test, n * n, prec);

        if (!_acb_vec_overlaps(th2, test, n * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 5);
            flint_printf("input:\n");
            _acb_vec_printd(thz, n, 5);
            _acb_vec_printd(th0, n, 5);
            flint_printf("output:\n");
            _acb_vec_printd(th2, n * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th0, n);
        _acb_vec_clear(thz, n);
        _acb_vec_clear(th2, n * n);
        _acb_vec_clear(test, n * n);
        _arb_vec_clear(d0, n);
        _arb_vec_clear(d, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
