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

    flint_printf("ql_new_roots_aux....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: does not fail for nice tau, roots are stored correctly */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong nb_z = 1 + n_randint(state, 3);
        slong prec = 100;
        slong highprec = 1000;
        slong nb_steps = 1 + n_randint(state, 10);
        acb_mat_t tau;        
        acb_ptr r, z, t, th, x;
        arb_ptr dist;
        slong k, j, a;

        acb_mat_init(tau, g, g);
        r = _acb_vec_init(2 * nb_z * n * nb_steps);
        z = _acb_vec_init(g * nb_z);
        t = _acb_vec_init(g);
        th = _acb_vec_init(n);
        x = _acb_vec_init(g);
        dist = _arb_vec_init(nb_z * n);

        acb_siegel_randtest_nice(tau, state, highprec);
        for (k = g; k < g * nb_z; k++)
        {
            acb_urandom(&z[k], state, highprec); /* z starts with 0 */
        }

        for (k = 0; k < nb_z; k++)
        {
            acb_theta_ql_sqr_dists_a(dist + k * n, z + k * g, tau, prec);
        }

        acb_theta_ql_new_roots_aux(r, t, z, nb_z, dist, tau, nb_steps, prec, highprec);

        if (!acb_is_finite(&r[0]))
        {
            flint_printf("FAIL (indeterminate)\n");
            flint_printf("g = %wd, nb_z = %wd, nb_steps = %wd, prec = %wd, tau:\n",
                g, nb_z, nb_steps, prec);
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        k = n_randint(state, nb_steps);
        j = n_randint(state, nb_z);

        /* Test: roots for 2^k (j-th vector z + t) */
        acb_mat_scalar_mul_2exp_si(tau, tau, k);
        _acb_vec_add(x, z + j * g, t, g, highprec);
        _acb_vec_scalar_mul_2exp_si(x, x, g, k);
        
        for (a = 0; a < n; a++)
        {
            acb_theta_naive_ind(th + a, a << g, x, 1, tau, highprec);
        }

        if (!_acb_vec_overlaps(th, r + 2 * j * nb_steps * n + k * n, n))
        {
            flint_printf("FAIL (values)\n");
            flint_printf("g = %wd, nb_z = %wd, nb_steps = %wd, prec = %wd, tau:\n",
                g, nb_z, nb_steps, prec);
            acb_mat_printd(tau, 5);
            flint_printf("Values:\n");
            _acb_vec_printd(th, n, 10);
            flint_printf("\n");
            _acb_vec_printd(r + 2 * j * nb_steps * n + k * n, n, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(r, 2 * nb_z * n * nb_steps);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(t, g);
        _acb_vec_clear(th, n);
        _acb_vec_clear(x, g);
        _arb_vec_clear(dist, n * nb_z);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
