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

    flint_printf("dist_a0....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: find zero value when z = tau a/2 + real stuff */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        slong n = 1 << g;
        slong prec = ACB_THETA_LOW_PREC;
        slong hprec = 200;
        slong bits = n_randint(state, 5);
        acb_mat_t tau;
        acb_ptr z;
        arb_ptr dist;
        arb_t c;
        ulong a = n_randint(state, n);
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dist = _arb_vec_init(n);
        arb_init(c);

        acb_siegel_randtest_reduced(tau, state, hprec, bits);
        acb_theta_char_get_acb(z, a, g);
        acb_mat_vector_mul_col(z, tau, z, prec);
        for (k = 0; k < g; k++)
        {
            arb_urandom(c, state, prec);
            acb_add_arb(&z[k], &z[k], c, prec);
        }

        acb_theta_dist_a0(dist, z, tau, prec);

        if (!arb_contains_zero(&dist[a]))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, a = %wd, tau:\n", g, a);
            acb_mat_printd(tau, 5);
            flint_printf("distances:\n");
            _arb_vec_printn(dist, n, 5, 0);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _arb_vec_clear(dist, n);
        arb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
