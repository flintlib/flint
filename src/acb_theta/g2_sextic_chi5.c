/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_g2_sextic_chi5(acb_poly_t res, acb_t chi5, const acb_mat_t tau, slong prec)
{
    slong g = 2;
    slong n2 = 1 << (2 * g);
    slong nb = acb_theta_jet_nb(1, g);
    fmpz_mat_t mat;
    acb_mat_t w, c, cinv;
    acb_ptr zero, dth, th;
    acb_t det;
    slong k;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(w, g, g);
    acb_mat_init(c, g, g);
    acb_mat_init(cinv, g, g);
    dth = _acb_vec_init(n2 * nb);
    th = _acb_vec_init(n2);
    zero = _acb_vec_init(g);
    acb_init(det);

    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform_cocycle_inv(w, c, cinv, mat, tau, prec);

    if (acb_siegel_is_reduced(w, -10, prec))
    {
        acb_theta_jet_all_notransform(dth, zero, 1, w, 1, prec);

        for (k = 0; k < n2; k++)
        {
            acb_set(&th[k], &dth[k * nb]);
        }
        acb_theta_g2_chi3_6(res, dth, prec);
        acb_theta_g2_chi5(chi5, th, prec);
        acb_poly_scalar_div(res, res, chi5, prec);

        acb_theta_g2_detk_symj(res, cinv, res, -2, 6, prec);
        acb_mat_det(det, cinv, prec);
        acb_pow_ui(det, det, 5, prec);

        if (acb_theta_g2_character(mat) == 1)
        {
            acb_neg(det, det);
        }
        acb_mul(chi5, chi5, det, prec);
    }
    else
    {
        /* Use sum_bound to avoid returning NaN */
        arb_t c, rho;
        arb_init(c);
        arb_init(rho);

        acb_theta_sum_bound(c, rho, zero, tau, 0);
        for (k = 0; k < nb * n2; k++)
        {
            arb_zero_pm_one(acb_realref(&dth[k]));
            arb_zero_pm_one(acb_imagref(&dth[k]));
            acb_mul_arb(&dth[k], &dth[k], c, prec);
            if (k % nb != 0) /* order 1 */
            {
                acb_div_arb(&dth[k], &dth[k], rho, prec);
            }
        }

        for (k = 0; k < n2; k++)
        {
            acb_set(&th[k], &dth[k * nb]);
        }
        acb_theta_g2_chi3_6(res, dth, prec);
        acb_theta_g2_chi5(chi5, th, prec);

        arb_clear(c);
        arb_clear(rho);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(w);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    _acb_vec_clear(dth, n2 * nb);
    _acb_vec_clear(th, n2);
    _acb_vec_clear(zero, g);
    acb_clear(det);
}
