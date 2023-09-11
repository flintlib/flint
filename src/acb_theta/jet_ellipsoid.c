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
    slong lp = ACB_THETA_LOW_PREC;
    arf_t R2, eps;
    arb_t t;
    arb_mat_t cho, invt;
    arb_ptr offset;

    arf_init(R2);
    arf_init(eps);
    arb_mat_init(cho, g, g);
    arb_mat_init(invt, g, g);
    offset = _arb_vec_init(g);
    arb_init(t);

    acb_theta_eld_cho(cho, tau, prec);

    if (arb_mat_is_finite(cho))
    {
        /* Get offset and bound on leading factor */
        arb_mat_one(invt);
        arb_mat_solve_triu(invt, cho, invt, 0, prec);
        arb_mat_transpose(invt, invt);
        _acb_vec_get_imag(offset, z, g);
        arb_mat_vector_mul_col(offset, invt, offset, prec);
        arb_const_pi(t, prec);
        _arb_vec_scalar_mul(offset, offset, g, t, prec);
        arb_zero(u);
        arb_dot(u, u, 0, offset, 1, offset, 1, g, prec);
        arb_exp(u, u, prec);

        /* Get radius, fill ellipsoid */
        acb_theta_jet_naive_radius(R2, eps, offset, cho, ord, prec);
        acb_theta_eld_fill(E, cho, R2, offset, lp);
        arb_mul_arf(u, u, eps, prec);
    }
    else
    {
        /* Cannot compute cho, result will be nan */
        arb_mat_one(cho);
        arf_zero(R2);
        acb_theta_eld_fill(E, cho, R2, offset, lp);
        arb_indeterminate(u);
    }

    arf_clear(R2);
    arf_clear(eps);
    arb_mat_clear(cho);
    arb_mat_clear(invt);
    _arb_vec_clear(offset, g);
    arb_clear(t);
}
