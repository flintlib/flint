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

    flint_printf("agm_mul_tight....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: respects relative precision */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 500);
        slong bits = n_randint(state, 3);
        slong delta = 25;
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, th0, r;
        arb_ptr dist, dist0;
        arb_t x;
        arf_t m, eps;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        r = _acb_vec_init(n);
        th = _acb_vec_init(n);
        th0 = _acb_vec_init(n);
        dist = _arb_vec_init(n);
        dist0 = _arb_vec_init(n);
        arb_init(x);
        arf_init(m);
        arf_init(eps);

        /* Generate distances, not too crazy */
        acb_siegel_randtest_nice(tau, state, prec);
        acb_theta_dist_a0(dist0, z, tau, prec);
        for (k = 0; k < g; k++)
        {
            acb_randtest_precise(&z[k], state, prec, bits);
        }
        acb_theta_dist_a0(dist, z, tau, prec);

        /* Generate values */
        for (k = 0; k < n; k++)
        {
            arb_neg(x, &dist[k]);
            arb_exp(x, x, prec);
            acb_urandom(&th[k], state, prec);
            acb_mul_arb(&th[k], &th[k], x, prec);

            arb_neg(x, &dist0[k]);
            arb_exp(x, x, prec);
            acb_urandom(&th0[k], state, prec);
            acb_mul_arb(&th0[k], &th0[k], x, prec);
        }

        acb_theta_agm_mul_tight(r, th0, th, dist0, dist, g, prec);
        acb_theta_agm_rel_mag_err(m, eps, r, dist, n, prec);

        /* Test: m <= 1 and eps is not too small */
        if (arf_cmp_si(m, 1) > 0 || arf_cmp_2exp_si(eps, -prec + delta) > 0)
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
            acb_mat_printd(tau, 5);
            flint_printf("distances:\n");
            _arb_vec_printn(dist0, n, 5, 0);
            flint_printf("\n");
            _arb_vec_printn(dist, n, 5, 0);
            flint_printf("\n");
            flint_printf("values:\n");
            _acb_vec_printd(th0, n, 5);
            flint_printf("\n");
            _acb_vec_printd(th, n, 5);
            flint_printf("\n");
            flint_printf("result:\n");
            _acb_vec_printd(r, n, 5);
            flint_printf("\n");
            flint_printf("m, eps: ");
            arf_printd(m, 10);
            flint_printf(", ");
            arf_printd(eps, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(r, n);
        _acb_vec_clear(th, n);
        _acb_vec_clear(th0, n);
        _arb_vec_clear(dist, n);
        _arb_vec_clear(dist0, n);
        arb_clear(x);
        arf_clear(m);
        arf_clear(eps);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
