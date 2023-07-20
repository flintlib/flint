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

    flint_printf("ql_new_roots....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: does not fail for z = 0, g <= 2 and nice tau */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n = 1 << g;
        slong prec = 100;
        slong nb_steps = n_randint(state, 10);
        acb_mat_t tau;
        acb_ptr r, z;
        arb_ptr dist;
        int res;
        
        acb_mat_init(tau, g, g);            
        r = _acb_vec_init(n * nb_steps);
        z = _acb_vec_init(g);
        dist = _arb_vec_init(n);
        
        acb_siegel_randtest_nice(tau, state, prec);
        acb_theta_ql_sqr_dists_a(dist, z, tau, prec);
        res = acb_theta_ql_new_roots(r, z, dist, tau, nb_steps, prec);

        if (!res)
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 10);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(r, n * nb_steps);
        _acb_vec_clear(z, g);
        _arb_vec_clear(dist, n);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
