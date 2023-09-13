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

    flint_printf("g2_covariants....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with g2_psi4 using psi4 = -(Co20 - 3*Co40)/20 */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 100);
        slong g = 2;
        slong n = 1 << (2 * g);
        acb_mat_t tau;
        acb_ptr z, th2;
        acb_poly_struct* r;
        acb_poly_t u, v;
        acb_t psi4, test;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        r = flint_malloc(26 * sizeof(acb_poly_struct));
        for (k = 0; k < 26; k++)
        {
            acb_poly_init(&r[k]);
        }
        acb_poly_init(u);
        acb_poly_init(v);
        acb_init(psi4);
        acb_init(test);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_theta_all(th2, z, tau, 1, prec);
        acb_theta_g2_psi4(psi4, th2, prec);

        acb_theta_g2_covariants(r, tau, prec);
        acb_poly_set_si(u, -3);
        acb_poly_mul(u, u, &r[8], prec);
        acb_poly_mul(v, &r[1], &r[1], prec);
        acb_poly_add(u, u, v, prec);
        acb_poly_get_coeff_acb(test, u, 0);
        acb_div_si(test, test, -20, prec);

        if (!acb_overlaps(psi4, test))
        {
            flint_printf("FAIL\n");
            acb_mat_printd(tau, 5);
            flint_printf("psi4, test:\n");
            acb_printd(psi4, 10);
            flint_printf("\n");
            acb_printd(test, 10);
            flint_printf("\nu:\n");
            acb_poly_printd(u, 5);
            flint_printf("\ncovariants:\n");
            for (k = 0; k < 26; k++)
            {
                acb_poly_printd(&r[k], 5);
                flint_printf("\n");
            }
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        for (k = 0; k < 26; k++)
        {
            acb_poly_clear(&r[k]);
        }
        flint_free(r);
        acb_poly_clear(u);
        acb_poly_clear(v);
        acb_clear(psi4);
        acb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
