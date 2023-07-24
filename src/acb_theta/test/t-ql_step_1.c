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

    flint_printf("ql_step_1....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive_ind */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong lp = ACB_THETA_LOW_PREC;
        slong prec = 100;
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr r, test, th, th0, roots;
        arb_ptr dist, dist0;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        r = _acb_vec_init(n);
        test = _acb_vec_init(n);
        th = _acb_vec_init(n);
        th0 = _acb_vec_init(n);
        roots = _acb_vec_init(n);
        dist = _arb_vec_init(n);
        dist0 = _arb_vec_init(n);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);

        /* Get input at zero */
        acb_theta_dist_a0(dist0, z, tau, lp);
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_ind(&th0[k], k << g, z, 1, tau, prec);
        }

        /* Get input at z */
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_dist_a0(dist, z, tau, lp);
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_ind(&th[k], k << g, z, 1, tau, prec);
        }

        /* Get output at tau/2, z/2 */
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        _acb_vec_scalar_mul_2exp_si(z, z, g, -1);
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_ind(&test[k], k << g, z, 1, tau, prec);
            acb_set_round(&roots[k], &test[k], lp);
        }

        acb_theta_ql_step_1(r, th, th0, roots, dist, dist0, g, prec);

        if (!acb_is_finite(&r[0]) || !_acb_vec_overlaps(r, test, n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 5);
            flint_printf("input:\n");
            _acb_vec_printd(th, n, 5);
            flint_printf("\n");
            _acb_vec_printd(th0, n, 5);
            flint_printf("\n");
            flint_printf("output:\n");
            _acb_vec_printd(r, n, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(r, n);
        _acb_vec_clear(test, n);
        _acb_vec_clear(th, n);
        _acb_vec_clear(th0, n);
        _acb_vec_clear(roots, n);
        _arb_vec_clear(dist, n);
        _arb_vec_clear(dist0, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
