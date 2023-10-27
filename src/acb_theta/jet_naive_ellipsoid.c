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
acb_theta_jet_naive_ellipsoid(acb_theta_eld_t E, arb_t u, acb_srcptr z,
    const acb_mat_t tau, slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arf_t R2, eps;
    arb_mat_t C, Yinv;
    arb_ptr v;
    int b;

    arf_init(R2);
    arf_init(eps);
    arb_mat_init(C, g, g);
    arb_mat_init(Yinv, g, g);
    v = _arb_vec_init(g);

    acb_siegel_yinv(Yinv, tau, prec);
    acb_siegel_cho(C, tau, prec);

    /* Get offset and bound on leading factor */
    _acb_vec_get_imag(v, z, g);
    arb_mat_vector_mul_col(v, Yinv, v, prec);
    arb_mat_vector_mul_col(v, C, v, prec);
    arb_zero(u);
    arb_dot(u, u, 0, v, 1, v, 1, g, prec);
    arb_exp(u, u, prec);

    /* Get radius, fill ellipsoid */
    acb_theta_jet_naive_radius(R2, eps, C, v, ord, prec);
    b = acb_theta_eld_set(E, C, R2, v);
    if (!b)
    {
        arb_pos_inf(u);
    }
    arb_mul_arf(u, u, eps, prec);

    arf_clear(R2);
    arf_clear(eps);
    arb_mat_clear(C);
    arb_mat_clear(Yinv);
    _arb_vec_clear(v, g);
}
