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

    flint_printf("g2_slash_basic_covariants....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with g2_psi4 using psi4 = -(Co20 - 3*Co40)/20 */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 100);
        slong g = 2;
        slong nb = ACB_THETA_G2_BASIC_NB;
        fmpz_mat_t mat;
        acb_mat_t tau, w, c;
        acb_poly_struct *cov1, *cov2, *test;
        acb_poly_t r1, r2;
        slong k;

        fmpz_mat_init(mat, 4, 4);
        acb_mat_init(tau, g, g);
        acb_mat_init(w, g, g);
        acb_mat_init(c, g, g);
        cov1 = flint_malloc(nb * sizeof(acb_poly_struct));
        cov2 = flint_malloc(nb * sizeof(acb_poly_struct));
        test = flint_malloc(nb * sizeof(acb_poly_struct));
        for (k = 0; k < nb; k++)
        {
            acb_poly_init(&cov1[k]);
            acb_poly_init(&cov2[k]);
            acb_poly_init(&test[k]);
        }
        acb_poly_init(r1);
        acb_poly_init(r2);

        acb_siegel_randtest_nice(tau, state, prec);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);
        acb_siegel_reduce(w, mat, tau, prec);
        acb_siegel_cocycle(c, mat, tau, prec);

        acb_theta_g2_fundamental_covariant(r1, tau, prec);
        acb_theta_g2_basic_covariants(cov1, r1, prec);
        acb_theta_g2_fundamental_covariant(r2, w, prec);
        acb_theta_g2_basic_covariants(cov2, r2, prec);
        acb_theta_g2_slash_basic_covariants(test, c, cov1, prec);

        for (k = 0; k < nb; k++)
        {
            if (!acb_poly_overlaps(&test[k], &cov2[k]))
            {
                flint_printf("FAIL (k = %wd)\n", k);
                fmpz_mat_print_pretty(mat);
                flint_printf("\n");
                acb_mat_printd(c, 5);
                acb_poly_printd(&test[k], 5);
                flint_printf("\n");
                acb_poly_printd(&cov2[k], 5);
                flint_printf("\nsource:\n");
                acb_poly_printd(&cov1[k], 5);
                flint_printf("\nfundamental:\n");
                acb_poly_printd(r1, 5);
                flint_printf("\n");
                acb_poly_printd(r2, 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        fmpz_mat_clear(mat);
        acb_mat_clear(tau);
        acb_mat_clear(w);
        acb_mat_clear(c);
        for (k = 0; k < nb; k++)
        {
            acb_poly_clear(&cov1[k]);
            acb_poly_clear(&cov2[k]);
            acb_poly_clear(&test[k]);
        }
        flint_free(cov1);
        flint_free(cov2);
        flint_free(test);
        acb_poly_clear(r1);
        acb_poly_clear(r2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
