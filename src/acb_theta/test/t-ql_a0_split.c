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

    flint_printf("ql_a0_split....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with ql_a0_naive using ql_a0_naive as worker */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 3);
        slong n = 1 << g;
        slong d = 1 + n_randint(state, g - 1);
        int has_t = iter % 2;
        slong nb_z = (has_t ? 3 : 1);
        slong prec = 50 + n_randint(state, 50);
        slong hprec = prec + 25;
        slong guard = 0;
        slong lp = ACB_THETA_LOW_PREC;
        acb_mat_t tau;
        acb_ptr z, t, r, test;
        arb_ptr dist, dist0;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        r = _acb_vec_init(nb_z * n);
        test = _acb_vec_init(2 * nb_z * n);
        dist = _arb_vec_init(n);
        dist0 = _arb_vec_init(n);

        acb_siegel_randtest_nice(tau, state, hprec);
        acb_theta_dist_a0(dist, z, tau, lp);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, hprec);
            if (has_t)
            {
                arb_urandom(acb_realref(&t[k]), state, hprec);
            }
        }
        acb_theta_dist_a0(dist, z, tau, lp);

        res = acb_theta_ql_a0_split(r, t, z, dist, tau, d, guard, prec,
            &acb_theta_ql_a0_naive);
        acb_theta_ql_a0_naive(test, t, z, dist0, dist, tau, guard, hprec);

        if (!_acb_vec_is_zero(z, g))
        {
            _acb_vec_set(test, test + nb_z * n, nb_z * n);
        }

        if (res && !_acb_vec_overlaps(r, test, nb_z * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 5);
            flint_printf("output:\n");
            _acb_vec_printd(r, nb_z * n, 5);
            _acb_vec_printd(test, nb_z * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r, nb_z * n);
        _acb_vec_clear(test, 2 * nb_z * n);
        _arb_vec_clear(dist, n);
        _arb_vec_clear(dist0, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
