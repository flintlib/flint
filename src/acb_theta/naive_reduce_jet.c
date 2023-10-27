/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_naive_reduce_jet(arb_ptr v, arb_t u, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t C, Yinv;

    arb_mat_init(C, g, g);
    arb_mat_init(Yinv, g, g);

    acb_siegel_cho(C, tau, prec);
    acb_siegel_yinv(Yinv, tau, prec);

    _acb_vec_get_imag(v, z, g);
    arb_mat_vector_mul_col(v, Yinv, v, prec);
    arb_mat_vector_mul_col(v, C, v, prec);
    arb_dot(u, NULL, 0, v, 1, v, 1, g, prec);
    arb_exp(u, u, prec);

    arb_mat_clear(C);
    arb_mat_clear(Yinv);
}
