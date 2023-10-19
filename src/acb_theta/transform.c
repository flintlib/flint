/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_transform_scal(acb_t scal, const fmpz_mat_t mat, acb_srcptr z,
    const acb_mat_t tau, slong kappa, slong prec)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t gamma;
    acb_mat_t w;
    acb_ptr Nz, v;
    acb_t x;

    fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
    acb_mat_init(w, g, g);
    v = _acb_vec_init(g);
    Nz = _acb_vec_init(g);
    acb_init(x);

    acb_one(x);
    acb_mul_2exp_si(x, x, -2);
    acb_exp_pi_i(x, x, prec);
    acb_pow_si(scal, x, kappa, prec);

    acb_siegel_transform_z(Nz, w, mat, z, tau, prec);
    acb_mat_set_fmpz_mat(w, gamma);
    acb_mat_vector_mul_col(v, w, z, prec);

    acb_dot(x, NULL, 0, v, 1, Nz, 1, g, prec);
    acb_exp_pi_i(x, x, prec);
    acb_mul(scal, scal, x, prec);

    fmpz_mat_window_clear(gamma);
    acb_mat_clear(w);
    _acb_vec_clear(v, g);
    _acb_vec_clear(Nz, g);
    acb_clear(x);
}

void
acb_theta_transform(acb_ptr res, const fmpz_mat_t mat, acb_srcptr th, acb_srcptr z,
    const acb_mat_t tau, int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong kappa;
    acb_mat_t c;
    acb_t scal, x;

    acb_mat_init(c, g, g);
    acb_init(scal);
    acb_init(x);

    if (sqr)
    {
        kappa = acb_theta_transform_kappa2(mat);
        acb_siegel_cocycle(c, mat, tau, prec);
        acb_mat_det(x, c, prec);
    }
    else
    {
        kappa = acb_theta_transform_kappa_new(x, mat, tau, prec);
    }
    acb_theta_transform_scal(scal, mat, z, tau, kappa, prec);
    if (sqr)
    {
        acb_sqr(scal, scal, prec);
    }
    acb_mul(scal, scal, x, prec);

    acb_theta_transform_proj(res, mat, th, sqr, prec);
    _acb_vec_scalar_mul(res, res, n * n, scal, prec);

    acb_mat_clear(c);
    acb_clear(scal);
    acb_clear(x);
}
