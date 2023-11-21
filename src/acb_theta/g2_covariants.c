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
acb_theta_g2_transvectants(acb_poly_struct * res, const acb_poly_t f, slong prec)
{
    acb_poly_t s;

    acb_poly_init(s);

    acb_poly_set(&res[0], f);
    acb_theta_g2_transvectant(&res[1], f, f, 6, 6, 6, prec);
    acb_theta_g2_transvectant(&res[2], f, f, 6, 6, 4, prec);
    acb_theta_g2_transvectant(&res[3], f, f, 6, 6, 2, prec);
    acb_theta_g2_transvectant(&res[4], f, &res[2], 6, 4, 4, prec);
    acb_theta_g2_transvectant(&res[5], f, &res[2], 6, 4, 2, prec);
    acb_theta_g2_transvectant(&res[6], f, &res[2], 6, 4, 1, prec);
    acb_theta_g2_transvectant(&res[7], f, &res[3], 6, 8, 1, prec);
    acb_theta_g2_transvectant(&res[8], &res[2], &res[2], 4, 4, 4, prec);
    acb_theta_g2_transvectant(&res[9], f, &res[4], 6, 2, 2, prec);
    acb_theta_g2_transvectant(&res[10], f, &res[4], 6, 2, 1, prec);
    acb_theta_g2_transvectant(&res[11], &res[3], &res[2], 8, 4, 1, prec);
    acb_theta_g2_transvectant(&res[12], &res[2], &res[4], 4, 2, 2, prec);
    acb_theta_g2_transvectant(&res[13], &res[2], &res[4], 4, 2, 1, prec);
    acb_theta_g2_transvectant(&res[14], &res[3], &res[4], 8, 2, 1, prec);
    acb_theta_g2_transvectant(&res[15], &res[4], &res[4], 2, 2, 2, prec);
    acb_theta_g2_transvectant(&res[16], &res[5], &res[4], 6, 2, 1, prec);
    acb_theta_g2_transvectant(&res[17], &res[6], &res[4], 8, 2, 2, prec);
    acb_poly_mul(s, &res[4], &res[4], prec); /* now C_32^2 */
    acb_theta_g2_transvectant(&res[18], f, s, 6, 4, 4, prec);
    acb_theta_g2_transvectant(&res[19], f, s, 6, 4, 3, prec);
    acb_theta_g2_transvectant(&res[20], &res[2], s, 4, 4, 3, prec);
    acb_theta_g2_transvectant(&res[21], &res[6], s, 8, 4, 4, prec);
    acb_poly_mul(s, s, &res[4], prec); /* now C_32^3 */
    acb_theta_g2_transvectant(&res[22], f, s, 6, 6, 6, prec);
    acb_theta_g2_transvectant(&res[23], f, s, 6, 6, 5, prec);
    acb_theta_g2_transvectant(&res[24], &res[6], s, 8, 6, 6, prec);
    acb_poly_mul(s, s, &res[4], prec); /* now C_32^4 */
    acb_theta_g2_transvectant(&res[25], &res[6], s, 8, 8, 8, prec);

    acb_poly_clear(s);
}

void
acb_theta_g2_covariants(acb_poly_struct * res, const acb_poly_t f, slong prec)
{
    double cofactors[ACB_THETA_G2_COV_NB] = {1, 60, 75, 90, 2250, 2250, 450,
        540, 11250, 67500, 13500, 13500, 168750, 67500, 405000, 10125000,
        2025000, 2700000, 151875000, 60750000, 15187500, 9112500000,
        227812500000, 13668750000, 8201250000000, 384433593750};
    acb_t c;
    fmpz_t m;
    slong k;

    acb_init(c);
    fmpz_init(m);

    acb_theta_g2_transvectants(res, f, prec);
    for (k = 0; k < ACB_THETA_G2_COV_NB; k++)
    {
        fmpz_set_d(m, cofactors[k]);
        acb_set_fmpz(c, m);
        acb_poly_scalar_mul(&res[k], &res[k], c, prec);
    }

    acb_clear(c);
    fmpz_clear(m);
}
