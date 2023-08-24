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

    flint_printf("siegel_cocycle_sqrtdet....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 6);
        fmpz_mat_t m1, m2, m3;
        slong k1, k2, k3;
        acb_mat_t tau1, tau2;
        acb_t c1, c2, c3, t;
        slong prec = 100 + n_randint(state, 200);
        slong mag_bits = n_randint(state, 5);

        fmpz_mat_init(m1, 2 * g, 2 * g);
        fmpz_mat_init(m2, 2 * g, 2 * g);
        fmpz_mat_init(m3, 2 * g, 2 * g);
        acb_mat_init(tau1, g, g);
        acb_mat_init(tau2, g, g);
        acb_init(c1);
        acb_init(c2);
        acb_init(c3);
        acb_init(t);

        acb_siegel_randtest(tau1, state, prec, mag_bits);
        acb_siegel_randtest(tau2, state, prec, mag_bits);
        sp2gz_randtest(m1, state, mag_bits);
        sp2gz_randtest(m2, state, mag_bits);
        fmpz_mat_mul(m3, m2, m1);

        /*k1 = acb_theta_transform_kappa(m1);
        k2 = acb_theta_transform_kappa(m2);
        k3 = acb_theta_transform_kappa(m3);*/
        k1 = 0; k2 = 0; k3 = 0;

        /* Test: chain rule */
        acb_siegel_cocycle_sqrtdet(c1, m1, tau1, prec);
        acb_siegel_transform(tau2, m1, tau1, prec);
        acb_siegel_cocycle_sqrtdet(c2, m2, tau2, prec);
        acb_siegel_cocycle_sqrtdet(c3, m3, tau1, prec);
        acb_mul(t, c2, c1, prec);

        if ((k1 + k2) % 8 != k3)
        {
            acb_neg(t, t);
        }

        if (!acb_overlaps(t, c3))
        {
            flint_printf("FAIL\n");
            acb_printd(c3, 10);
            flint_printf("\n");
            acb_printd(t, 10);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(m1);
        fmpz_mat_clear(m2);
        fmpz_mat_clear(m3);
        acb_mat_clear(tau1);
        acb_mat_clear(tau2);
        acb_clear(c1);
        acb_clear(c2);
        acb_clear(c3);
        acb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
