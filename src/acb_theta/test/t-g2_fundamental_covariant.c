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

    flint_printf("g2_fundamental_covariant....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with chi6m2 up to Pi^6 */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong n = 1 << (2 * g);
        slong nb = acb_theta_jet_nb(1, g + 1);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_mat_t tau;
        acb_ptr z, dth;
        acb_poly_t chi, test;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        dth = _acb_vec_init(n * nb);
        acb_poly_init(chi);
        acb_poly_init(test);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);

        acb_theta_jet_all(dth, z, tau, 1, prec);
        acb_theta_g2_chi6m2(test, dth, prec);
        acb_theta_g2_fundamental_covariant(chi, tau, prec);

        if (!acb_poly_overlaps(chi, test))
        {
            flint_printf("FAIL\n");
            flint_printf("chi:\n");
            acb_poly_printd(chi, 5);
            flint_printf("\ntest:\n");
            acb_poly_printd(test, 5);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(dth, n * nb);
        acb_poly_clear(chi);
        acb_poly_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
