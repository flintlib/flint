/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_theta.h"

static void
acb_theta_g2_transvectants(acb_ptr res, const acb_poly_t f, slong prec)
{
    acb_poly_t s, r2, r3, r4, r5, r6;

    acb_poly_init(s);
    acb_poly_init(r2);
    acb_poly_init(r3);
    acb_poly_init(r4);
    acb_poly_init(r5);
    acb_poly_init(r6);

    /* Get polynomials */
    acb_theta_g2_transvectant(r2, f, f, 6, 6, 4, prec);
    acb_theta_g2_transvectant(r3, f, f, 6, 6, 2, prec);
    acb_theta_g2_transvectant(r4, f, r2, 6, 4, 4, prec);
    acb_theta_g2_transvectant(r5, f, r2, 6, 4, 2, prec);
    acb_theta_g2_transvectant(r6, f, r2, 6, 4, 1, prec);

    /* Get leading coefficients of f, r2, ..., r6 */
    acb_poly_get_coeff_acb(&res[0], f, 6);
    acb_poly_get_coeff_acb(&res[2], r2, 4);
    acb_poly_get_coeff_acb(&res[3], r3, 8);
    acb_poly_get_coeff_acb(&res[4], r4, 2);
    acb_poly_get_coeff_acb(&res[5], r5, 6);
    acb_poly_get_coeff_acb(&res[6], r6, 8);

    /* Get other coefficients */
    acb_theta_g2_transvectant_lead(&res[1], f, f, 6, 6, 6, prec);
    acb_theta_g2_transvectant_lead(&res[7], f, r3, 6, 8, 1, prec);
    acb_theta_g2_transvectant_lead(&res[8], r2, r2, 4, 4, 4, prec);
    acb_theta_g2_transvectant_lead(&res[9], f, r4, 6, 2, 2, prec);
    acb_theta_g2_transvectant_lead(&res[10], f, r4, 6, 2, 1, prec);
    acb_theta_g2_transvectant_lead(&res[11], r3, r2, 8, 4, 1, prec);
    acb_theta_g2_transvectant_lead(&res[12], r2, r4, 4, 2, 2, prec);
    acb_theta_g2_transvectant_lead(&res[13], r2, r4, 4, 2, 1, prec);
    acb_theta_g2_transvectant_lead(&res[14], r3, r4, 8, 2, 1, prec);
    acb_theta_g2_transvectant_lead(&res[15], r4, r4, 2, 2, 2, prec);
    acb_theta_g2_transvectant_lead(&res[16], r5, r4, 6, 2, 1, prec);
    acb_theta_g2_transvectant_lead(&res[17], r6, r4, 8, 2, 2, prec);
    acb_poly_mul(s, r4, r4, prec); /* C_32^2 */
    acb_theta_g2_transvectant_lead(&res[18], f, s, 6, 4, 4, prec);
    acb_theta_g2_transvectant_lead(&res[19], f, s, 6, 4, 3, prec);
    acb_theta_g2_transvectant_lead(&res[20], r2, s, 4, 4, 3, prec);
    acb_theta_g2_transvectant_lead(&res[21], r6, s, 8, 4, 4, prec);
    acb_poly_mul(s, s, r4, prec); /* now C_32^3 */
    acb_theta_g2_transvectant_lead(&res[22], f, s, 6, 6, 6, prec);
    acb_theta_g2_transvectant_lead(&res[23], f, s, 6, 6, 5, prec);
    acb_theta_g2_transvectant_lead(&res[24], r6, s, 8, 6, 6, prec);
    acb_poly_mul(s, s, r4, prec); /* now C_32^4 */
    acb_theta_g2_transvectant_lead(&res[25], r6, s, 8, 8, 8, prec);

    acb_poly_clear(s);
    acb_poly_clear(r2);
    acb_poly_clear(r3);
    acb_poly_clear(r4);
    acb_poly_clear(r5);
    acb_poly_clear(r6);
}

void
acb_theta_g2_covariants_lead(acb_ptr res, const acb_poly_t f, slong prec)
{
    double cofactors[ACB_THETA_G2_COV_NB] = {1, 60, 75, 90, 2250, 2250, 450,
        540, 11250, 67500, 13500, 13500, 168750, 67500, 405000, 10125000,
        2025000, 2700000, 151875000, 60750000, 15187500, 9112500000,
        227812500000, 13668750000, 8201250000000, 384433593750};
    fmpz_t m;
    slong k;

    fmpz_init(m);

    acb_theta_g2_transvectants(res, f, prec);
    for (k = 0; k < ACB_THETA_G2_COV_NB; k++)
    {
        fmpz_set_d(m, cofactors[k]);
        acb_mul_fmpz(&res[k], &res[k], m, prec);
    }

    fmpz_clear(m);
}
