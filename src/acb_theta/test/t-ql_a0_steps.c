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

    flint_printf("ql_a0_steps....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with ql_a0_naive using ql_a0_naive as worker */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 2);
        slong n = 1 << g;
        slong d = 1 + n_randint(state, g);
        int has_t = iter % 2;
        int has_z = (iter % 4) / 2;
        slong nbt = (has_t ? 3 : 1);
        slong nbz = (has_z ? 2 : 1);
        slong prec = 200 + n_randint(state, 1000);
        slong hprec = prec + 50;
        slong guard = ACB_THETA_LOW_PREC;
        slong lp = ACB_THETA_LOW_PREC;
        acb_mat_t tau;
        acb_ptr z, zero, t, r, test;
        arb_ptr dist, dist0;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        zero = _acb_vec_init(g);
        t = _acb_vec_init(g);
        r = _acb_vec_init(nbz * nbt * n);
        test = _acb_vec_init(nbz * nbt * n);
        dist = _arb_vec_init(n);
        dist0 = _arb_vec_init(n);

        acb_siegel_randtest_nice(tau, state, hprec);
        for (k = d; k < g; k++)
        {
            for (j = d; j < g; j++)
            {
                acb_mul_2exp_si(acb_mat_entry(tau, j, k),
                    acb_mat_entry(tau, j, k), 2 * ACB_THETA_QL_SPLIT + 1);
            }
        }
        for (k = 0; k < g; k++)
        {
            if (has_z)
            {
                acb_urandom(&z[k], state, hprec);
            }
            if (has_t)
            {
                arb_urandom(acb_realref(&t[k]), state, hprec);
            }
        }
        acb_theta_dist_a0(dist, z, tau, lp);
        acb_theta_dist_a0(dist0, zero, tau, lp);

        res = acb_theta_ql_a0_steps(r, t, z, dist0, dist, tau, guard, prec,
            &acb_theta_ql_a0_naive);
        acb_theta_ql_a0_naive(test, t, z, dist0, dist, tau, guard, hprec);

        if (res && !_acb_vec_overlaps(r, test, nbz * nbt * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, d = %wd, has_z = %wd, has_t = %wd, tau:\n",
                g, prec, d, has_z, has_t);
            acb_mat_printd(tau, 5);
            flint_printf("output:\n");
            _acb_vec_printd(r, nbz * nbt * n, 5);
            flint_printf("\n");
            _acb_vec_printd(test, nbz * nbt * n, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(zero, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r, nbz * nbt * n);
        _acb_vec_clear(test, nbz * nbt * n);
        _arb_vec_clear(dist, n);
        _arb_vec_clear(dist0, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
