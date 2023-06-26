/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_transform_scal(acb_t scal_z, acb_t scal_0, acb_srcptr z,
                         const acb_mat_t tau, const fmpz_mat_t mat, slong k2,
                         slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t c;
    acb_mat_t w;
    acb_mat_t vec;
    acb_ptr Nz;
    acb_t mu;
    acb_t det;
    slong k;

    fmpz_mat_init(c, g, g);
    acb_mat_init(w, g, g);
    acb_mat_init(vec, g, 1);
    Nz = _acb_vec_init(g);
    acb_init(mu);
    acb_init(det);

    acb_onei(mu);
    acb_pow_si(mu, mu, k2, prec);
    acb_siegel_cocycle(w, mat, tau, prec);
    acb_mat_det(det, w, prec);
    acb_mul(scal_0, det, mu, prec);

    acb_siegel_transform_ext(Nz, w, mat, z, tau, prec);
    fmpz_mat_get_c(c, mat);
    acb_mat_set_fmpz_mat(w, c);
    for (k = 0; k < g; k++)
    {
        acb_set(acb_mat_entry(vec, k, 0), &z[k]);
    }
    acb_mat_mul(vec, w, vec, prec);
    acb_zero(det);
    for (k = 0; k < g; k++)
    {
        acb_addmul(det, &Nz[k], acb_mat_entry(vec, k, 0), prec);
    }
    acb_mul_2exp_si(det, det, 1);
    acb_exp_pi_i(det, det, prec);
    acb_mul(scal_z, scal_0, det, prec);

    fmpz_mat_clear(c);
    acb_mat_clear(w);
    acb_mat_clear(vec);
    _acb_vec_clear(Nz, g);
    acb_clear(mu);
    acb_clear(det);
}
