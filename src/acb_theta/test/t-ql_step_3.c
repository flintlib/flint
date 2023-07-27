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

    flint_printf("ql_step_3....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive_ind */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong lp = ACB_THETA_LOW_PREC;
        slong prec = 100;
        acb_mat_t tau;
        acb_ptr z, t, x;
        acb_ptr r, test, th, th0, roots;
        arb_ptr dist, dist0;
        slong j, k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        x = _acb_vec_init(g);
        r = _acb_vec_init(3 * n);
        test = _acb_vec_init(3 * n);
        th = _acb_vec_init(3 * n);
        th0 = _acb_vec_init(3 * n);
        roots = _acb_vec_init(2 * n);
        dist = _arb_vec_init(n);
        dist0 = _arb_vec_init(n);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, 1);
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec);
        }

        /* Get input at zero */
        acb_theta_dist_a0(dist0, z, tau, lp);
        for (j = 0; j < 3; j++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, j, prec);
            for (k = 0; k < n; k++)
            {
                acb_theta_naive_ind(&th0[j * n + k], k << g, x, 1, tau, prec);
            }
        }

        /* Get input at z */
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }
        acb_theta_dist_a0(dist, z, tau, lp);
        for (j = 0; j < 3; j++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, j, prec);
            _acb_vec_add(x, x, z, g, prec);
            for (k = 0; k < n; k++)
            {
                acb_theta_naive_ind(&th[j * n + k], k << g, x, 1, tau, prec);
            }
        }

        /* Get output at tau/2, z/2 */
        acb_mat_scalar_mul_2exp_si(tau, tau, -1);
        _acb_vec_scalar_mul_2exp_si(z, z, g, -1);
        _acb_vec_scalar_mul_2exp_si(t, t, g, -1);
        for (j = 0; j < 3; j++)
        {
            _acb_vec_scalar_mul_ui(x, t, g, j, prec);
            _acb_vec_add(x, x, z, g, prec);
            for (k = 0; k < n; k++)
            {
                acb_theta_naive_ind(&test[j * n + k], k << g, x, 1, tau, prec);
                if (j > 0)
                {
                    acb_set_round(&roots[(j - 1) * n + k], &test[j * n + k], lp);
                }
            }
        }

        acb_theta_ql_step_3(r, th0, th, roots, dist0, dist, g, prec);

        if (!acb_is_finite(&r[0]) || !_acb_vec_overlaps(r, test, 3 * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 5);
            flint_printf("input:\n");
            _acb_vec_printd(th, 3 * n, 5);
            _acb_vec_printd(th0, 3 * n, 5);
            flint_printf("output:\n");
            _acb_vec_printd(r, 3 * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(x, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r, 3 * n);
        _acb_vec_clear(test, 3 * n);
        _acb_vec_clear(th, 3 * n);
        _acb_vec_clear(th0, 3 * n);
        _acb_vec_clear(roots, 2 * n);
        _arb_vec_clear(dist, n);
        _arb_vec_clear(dist0, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
