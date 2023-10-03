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
    arb_mat_t cho, Yinv;
    arb_ptr offset;

    arf_init(R2);
    arf_init(eps);
    arb_mat_init(cho, g, g);
    arb_mat_init(Yinv, g, g);
    offset = _arb_vec_init(g);

    acb_theta_eld_cho(cho, tau, prec);
    acb_mat_get_imag(Yinv, tau);
    arb_mat_inv(Yinv, Yinv, prec);

    if (arb_mat_is_finite(cho))
    {
        /* Get offset and bound on leading factor */
        _acb_vec_get_imag(offset, z, g);
        arb_mat_vector_mul_col(offset, Yinv, offset, prec);
        arb_mat_vector_mul_col(offset, cho, offset, prec);
        arb_zero(u);
        arb_dot(u, u, 0, offset, 1, offset, 1, g, prec);
        arb_exp(u, u, prec);

        /* Get radius, fill ellipsoid */
        acb_theta_jet_naive_radius(R2, eps, offset, cho, ord, prec);
        acb_theta_eld_fill(E, cho, R2, offset);
        arb_mul_arf(u, u, eps, prec);
    }
    else
    {
        /* Cannot compute cho, result will be nan */
        arb_mat_one(cho);
        arf_zero(R2);
        acb_theta_eld_fill(E, cho, R2, offset);
        arb_indeterminate(u);
    }

    arf_clear(R2);
    arf_clear(eps);
    arb_mat_clear(cho);
    arb_mat_clear(Yinv);
    _arb_vec_clear(offset, g);
}
