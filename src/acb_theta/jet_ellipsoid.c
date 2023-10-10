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
acb_theta_jet_ellipsoid(acb_theta_eld_t E, arb_t u, acb_srcptr z,
    const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arf_t R2, eps;
    arb_mat_t C, Yinv;
    arb_ptr v;

    arf_init(R2);
    arf_init(eps);
    arb_mat_init(C, g, g);
    arb_mat_init(Yinv, g, g);
    v = _arb_vec_init(g);

    acb_theta_eld_cho(C, tau, prec);
    acb_mat_get_imag(Yinv, tau);
    arb_mat_inv(Yinv, Yinv, prec);

    if (arb_mat_is_finite(C))
    {
        /* Get offset and bound on leading factor */
        _acb_vec_get_imag(v, z, g);
        arb_mat_vector_mul_col(v, Yinv, v, prec);
        arb_mat_vector_mul_col(v, C, v, prec);
        arb_zero(u);
        arb_dot(u, u, 0, v, 1, v, 1, g, prec);
        arb_exp(u, u, prec);

        /* Get radius, fill ellipsoid */
        acb_theta_jet_naive_radius(R2, eps, C, v, ord, prec);
        acb_theta_eld_fill(E, C, R2, v);
        arb_mul_arf(u, u, eps, prec);
    }
    else
    {
        /* Cannot compute C, result will be nan */
        arb_mat_one(C);
        arf_zero(R2);
        acb_theta_eld_fill(E, C, R2, v);
        arb_indeterminate(u);
    }

    arf_clear(R2);
    arf_clear(eps);
    arb_mat_clear(C);
    arb_mat_clear(Yinv);
    _arb_vec_clear(v, g);
}
