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

    flint_printf("transform_kappa_new....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: matches combination of transform_kappa and transform_sqrtdet */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong bits = 0;
        slong prec = 200;
        fmpz_mat_t mat;
        fmpz_mat_t x;
        acb_mat_t tau;
        acb_t sqrtdet, test;
        slong kappa1, kappa2;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        fmpz_mat_init(x, 2, 2);
        acb_mat_init(tau, g, g);
        acb_init(sqrtdet);
        acb_init(test);

        sp2gz_randtest(mat, state, bits);
        acb_siegel_randtest_nice(tau, state, prec);
        
        flint_printf("\n\ntau, mat:\n");
        acb_mat_printd(tau, 5);
        fmpz_mat_print_pretty(mat);
        flint_printf("\n");

        kappa1 = acb_theta_transform_kappa(mat);
        acb_theta_transform_sqrtdet(test, mat, tau, prec);

        flint_printf("found kappa = %wd, sqrtdet = ", kappa1);
        acb_printd(test, 5);
        flint_printf("\n");

        kappa2 = acb_theta_transform_kappa_new(sqrtdet, mat, tau, prec);
        flint_printf("found new kappa = %wd, sqrtdet = ", kappa2);
        acb_printd(sqrtdet, 5);
        flint_printf("\n");

        if (kappa1 % 4 != kappa2 % 4)
        {
            flint_printf("FAIL (kappa)\n");
            flint_abort();
        }

        if (kappa1 != kappa2)
        {
            acb_neg(test, test);
        }
        if (!acb_overlaps(test, sqrtdet))
        {
            flint_printf("FAIL (sqrtdet)\n");
            flint_abort();
        }

        fmpz_mat_clear(mat);
        fmpz_mat_clear(x);
        acb_mat_clear(tau);
        acb_clear(sqrtdet);
        acb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

