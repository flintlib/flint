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

    flint_printf("uql_roots....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: does not fail for nice tau, roots are stored correctly */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong nb_z = 1 + n_randint(state, 3);
        slong prec = 500;
        slong nb_steps = 1 + n_randint(state, 3);
        acb_mat_t tau;
        acb_ptr r, z, t, th, x;
        slong res;
        slong k, j;

        acb_mat_init(tau, g, g);
        r = _acb_vec_init(2 * (nb_z + 1) * n * nb_steps);
        z = _acb_vec_init(g * nb_z);
        t = _acb_vec_init(g);
        th = _acb_vec_init(n);
        x = _acb_vec_init(g);

        acb_siegel_randtest_nice(tau, state, prec);
        for (k = 0; k < g * nb_z; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        res = acb_theta_uql_roots(r, t, z, nb_z, tau, nb_steps, prec);

        if (res < 0)
        {
            flint_printf("FAIL (res)\n");
            flint_printf("g = %wd, nb_z = %wd, nb_steps = %wd, prec = %wd, tau:\n",
                g, nb_z, nb_steps, prec);
            acb_mat_printd(tau, 5);
            flint_abort();
        }

        k = n_randint(state, nb_steps);
        j = n_randint(state, nb_z);

        /* Test: roots for 2^k t */
        acb_mat_scalar_mul_2exp_si(tau, tau, k);
        _acb_vec_scalar_mul_2exp_si(x, t, g, k);
        acb_theta_naive_a0(th, x, 1, tau, prec);

        if (!_acb_vec_overlaps(th, r + 2 * (nb_z + 1) * n * k, n))
        {
            flint_printf("FAIL (values at 2^k t)\n");
            flint_printf("g = %wd, nb_z = %wd, nb_steps = %wd, prec = %wd, tau:\n",
                g, nb_z, nb_steps, prec);
            acb_mat_printd(tau, 5);
            flint_printf("Values:\n");
            _acb_vec_printd(th, n, 10);
            flint_printf("\n");
            _acb_vec_printd(r + 2 * (nb_z + 1) * n * k, n, 10);
            flint_printf("\n");
            flint_abort();
        }

        /* Test: roots for 2^k(j-th vector z + 2t) */
        _acb_vec_scalar_mul_2exp_si(t, t, g, 1);
        _acb_vec_add(x, z + j * g, t, g, prec);
        _acb_vec_scalar_mul_2exp_si(x, x, g, k);
        acb_theta_naive_a0(th, x, 1, tau, prec);
        j = 2 * (nb_z + 1) * n * k + 2 * (j + 1) * n + n;

        if (!_acb_vec_overlaps(th, r + j, n))
        {
            flint_printf("FAIL (values at 2^k(z + 2t))\n");
            flint_printf("g = %wd, nb_z = %wd, nb_steps = %wd, prec = %wd, tau:\n",
                g, nb_z, nb_steps, prec);
            acb_mat_printd(tau, 5);
            flint_printf("Values:\n");
            _acb_vec_printd(th, n, 10);
            flint_printf("\n");
            _acb_vec_printd(r + j, n, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(r, 2 * (nb_z + 1) * n * nb_steps);
        _acb_vec_clear(z, g * nb_z);
        _acb_vec_clear(t, g);
        _acb_vec_clear(th, n);
        _acb_vec_clear(x, g);        
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
