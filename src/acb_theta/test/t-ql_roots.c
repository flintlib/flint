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

    flint_printf("ql_roots....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: does not fail for z = 0, g <= 2 and nice tau */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 2);
        slong n = 1 << g;
        slong prec = 1000 + n_randint(state, 1000);
        slong nb_steps;
        slong nb_z = n_randint(state, 10);
        acb_mat_t tau;
        acb_ptr r, z;
        slong res;
        
        acb_mat_init(tau, g, g);
        acb_siegel_randtest_nice(tau, state, prec);
        nb_steps = acb_theta_ql_nb_steps(tau, prec);
            
        r = _acb_vec_init(n * nb_steps * nb_z);
        z = _acb_vec_init(nb_z * g);

        res = acb_theta_ql_roots(r, z, nb_z, tau, nb_steps, prec);

        if (res <= 0)
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 10);
            flint_printf("nb_z = %wd, prec = %wd\n", nb_z, prec);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(r, n * nb_steps * nb_z);
        _acb_vec_clear(z, nb_z * g);
    }
    
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
