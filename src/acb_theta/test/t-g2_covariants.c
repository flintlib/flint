/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

#define ACB_THETA_G2_COV_K {1,2,2,2,3,3,3,3,4,4,4,4,5,5,5,6,6,6,7,7,8,9,10,10,12,15}
#define ACB_THETA_G2_COV_J {6,0,4,8,2,6,8,12,0,4,6,10,2,4,8,0,6,6,2,4,2,4,0,2,2,0}

TEST_FUNCTION_START(acb_theta_g2_covariants, state)
{
    slong iter;

    /* Test:
       - agrees with g2_psi4 using psi4 = -(Co20 - 3*Co40)/20
       - covariants transform as det^(k-j/2)Sym^j
       - covariants take integral values on integral polynomials */
    for (iter = 0; iter < 5 * flint_test_multiplier(); iter++)
    {
        slong prec = 200 + n_randint(state, 200);
        slong g = 2;
        slong n = 1 << (2 * g);
        slong bits = 2;
        slong jlist[26] = ACB_THETA_G2_COV_J;
        slong klist[26] = ACB_THETA_G2_COV_K;
        fmpz_mat_t mat;
        acb_mat_t tau, w, c;
        acb_ptr z, th2;
        acb_poly_struct * cov1;
        acb_poly_struct * cov2;
        acb_poly_t u, v;
        fmpz_poly_t pol;
        acb_t psi4, test;
        slong k;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        acb_mat_init(tau, g, g);
        acb_mat_init(w, g, g);
        acb_mat_init(c, g, g);
        z = _acb_vec_init(g);
        th2 = _acb_vec_init(n);
        cov1 = flint_malloc(26 * sizeof(acb_poly_struct));
        cov2 = flint_malloc(26 * sizeof(acb_poly_struct));
        for (k = 0; k < 26; k++)
        {
            acb_poly_init(&cov1[k]);
            acb_poly_init(&cov2[k]);
        }
        acb_poly_init(u);
        acb_poly_init(v);
        fmpz_poly_init(pol);
        acb_init(psi4);
        acb_init(test);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        sp2gz_randtest(mat, state, bits);

        acb_theta_all(th2, z, tau, 1, prec);
        acb_theta_g2_psi4(psi4, th2, prec);
        acb_theta_g2_sextic(u, tau, prec);
        acb_theta_g2_covariants(cov1, u, prec);

        acb_siegel_transform(w, mat, tau, prec);
        acb_siegel_cocycle(c, mat, tau, prec);
        acb_theta_g2_sextic(u, w, prec);
        acb_theta_g2_covariants(cov2, u, prec);

        /* Test psi4 */
        acb_poly_set_si(u, -3);
        acb_poly_mul(u, u, &cov1[8], prec);
        acb_poly_mul(v, &cov1[1], &cov1[1], prec);
        acb_poly_add(u, u, v, prec);
        acb_poly_get_coeff_acb(test, u, 0);
        acb_div_si(test, test, -20, prec);

        if (!acb_overlaps(psi4, test))
        {
            flint_printf("FAIL (psi4)\n");
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
                acb_poly_printd(&cov1[k], 5);
                flint_printf("\n");
            }
            flint_abort();
        }

        /* Test transformation */
        for (k = 0; k < 26; k++)
        {
            acb_theta_g2_detk_symj(u, c, &cov1[k], klist[k] - jlist[k]/2, jlist[k], prec);
            if (!acb_poly_overlaps(u, &cov2[k]))
            {
                flint_printf("FAIL (transform, k = %wd)\n", k);
                acb_mat_printd(tau, 5);
                fmpz_mat_print_pretty(mat);
                flint_printf("\n");
                acb_poly_printd(u, 5);
                flint_printf("\n");
                acb_poly_printd(&cov2[k], 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        /* Test integrality */
        acb_poly_zero(u);
        for (k = 0; k <= 6; k++)
        {
            acb_poly_set_coeff_si(u, k, n_randint(state, 10));
        }
        acb_theta_g2_covariants(cov1, u, prec);
        for (k = 0; k < 26; k++)
        {
            if (!acb_poly_get_unique_fmpz_poly(pol, &cov1[k]))
            {
                flint_printf("FAIL (integrality, k = %wd)\n", k);
                acb_poly_printd(&cov1[k], 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        fmpz_mat_clear(mat);
        acb_mat_clear(tau);
        acb_mat_clear(w);
        acb_mat_clear(c);
        _acb_vec_clear(z, g);
        _acb_vec_clear(th2, n);
        for (k = 0; k < 26; k++)
        {
            acb_poly_clear(&cov1[k]);
            acb_poly_clear(&cov2[k]);
        }
        flint_free(cov1);
        flint_free(cov2);
        acb_poly_clear(u);
        acb_poly_clear(v);
        fmpz_poly_clear(pol);
        acb_clear(psi4);
        acb_clear(test);
    }

    TEST_FUNCTION_END(state);
}
