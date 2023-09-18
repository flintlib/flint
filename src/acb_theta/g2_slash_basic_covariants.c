/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_g2_slash_basic_covariants(acb_poly_struct* res, const acb_mat_t c,
    const acb_poly_struct* cov, slong prec)
{
    slong klist[] = ACB_THETA_G2_BASIC_K;
    slong jlist[] = ACB_THETA_G2_BASIC_J;
    slong nb = ACB_THETA_G2_BASIC_NB;
    slong nb_j = ACB_THETA_G2_MAX_J/2 + 1;
    acb_t det;
    acb_ptr det_pow;
    acb_ptr inv_pow;
    acb_poly_t x, y;
    acb_poly_struct* pow_x;
    acb_poly_struct* pow_y;
    acb_poly_struct** products;
    slong i, j, k, e;

    /* Init everything */
    acb_init(det);
    det_pow = _acb_vec_init(ACB_THETA_G2_MAX_K + 1);
    inv_pow = _acb_vec_init(ACB_THETA_G2_NEG_EXP + 1);
    acb_poly_init(x);
    acb_poly_init(y);
    pow_x = flint_malloc((ACB_THETA_G2_MAX_J + 1) * sizeof(acb_poly_struct));
    pow_y = flint_malloc((ACB_THETA_G2_MAX_J + 1) * sizeof(acb_poly_struct));
    for (k = 0; k < ACB_THETA_G2_MAX_J + 1; k++)
    {
        acb_poly_init(&pow_x[k]);
        acb_poly_init(&pow_y[k]);
    }
    products = flint_malloc(nb_j * sizeof(acb_poly_struct*));
    for (k = 0; k < nb_j; k++)
    {
        products[k] = flint_malloc((2 * k + 1) * sizeof(acb_poly_struct));
        for (j = 0; j < 2 * k + 1; j++)
        {
            acb_poly_init(&products[k][j]);
        }
    }

    /* Precompute products and powers of det */
    acb_mat_det(det, c, prec);
    _acb_vec_set_powers(det_pow, det, ACB_THETA_G2_MAX_K + 1, prec);
    acb_inv(det, det, prec);
    _acb_vec_set_powers(inv_pow, det, ACB_THETA_G2_NEG_EXP + 1, prec);
    acb_poly_set_coeff_acb(x, 0, acb_mat_entry(c, 1, 0));
    acb_poly_set_coeff_acb(x, 1, acb_mat_entry(c, 0, 0));
    acb_poly_set_coeff_acb(y, 0, acb_mat_entry(c, 1, 1));
    acb_poly_set_coeff_acb(y, 1, acb_mat_entry(c, 0, 1));
    acb_poly_one(&pow_x[0]);
    acb_poly_one(&pow_y[0]);
    for (k = 1; k < ACB_THETA_G2_MAX_J + 1; k++)
    {
        acb_poly_mul(&pow_x[k], &pow_x[k - 1], x, prec);
        acb_poly_mul(&pow_y[k], &pow_y[k - 1], y, prec);
    }
    for (k = 0; k < nb_j; k++)
    {
        for (j = 0; j < 2 * k + 1; j++)
        {
            acb_poly_mul(&products[k][j], &pow_x[j], &pow_y[2 * k - j], prec);
        }
    }

    /* Make substitutions and scalar products */
    for (i = 0; i < nb; i++)
    {
        acb_theta_g2_subst_covariant(&res[i], products[jlist[i]/2], &cov[i], jlist[i], prec);
        e = klist[i] - jlist[i]/2;
        if (e >= 0)
        {
            acb_poly_scalar_mul(&res[i], &res[i], &det_pow[e], prec);
        }
        else
        {
            acb_poly_scalar_mul(&res[i], &res[i], &inv_pow[-e], prec);
        }
    }

    /* Clear */
    acb_clear(det);
    _acb_vec_clear(det_pow, ACB_THETA_G2_MAX_K + 1);
    _acb_vec_clear(inv_pow, ACB_THETA_G2_NEG_EXP + 1);
    acb_poly_clear(x);
    acb_poly_clear(y);
    for (k = 0; k < ACB_THETA_G2_MAX_J + 1; k++)
    {
        acb_poly_clear(&pow_x[k]);
        acb_poly_clear(&pow_y[k]);
    }
    flint_free(pow_x);
    flint_free(pow_y);
    for (k = 0; k < nb_j; k++)
    {
        for (j = 0; j < 2 * k + 1; j++)
        {
            acb_poly_clear(&products[k][j]);
        }
        flint_free(products[k]);
    }
    flint_free(products);
}
