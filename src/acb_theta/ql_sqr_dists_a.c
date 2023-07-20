/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_ql_sqr_dists_a(arb_ptr dist, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_mat_t Yinv, cho;
    arb_ptr y, offset;
    arb_t pi;
    ulong a;

    arb_mat_init(Yinv, g, g);
    arb_mat_init(cho, g, g);
    y = _arb_vec_init(g);
    offset = _arb_vec_init(g);
    arb_init(pi);

    arb_const_pi(pi, prec);
    acb_mat_get_imag(Yinv, tau);
    arb_mat_scalar_mul_arb(cho, Yinv, pi, prec);
    arb_mat_cho(cho, cho, prec);
    arb_mat_transpose(cho, cho);
    arb_mat_inv(Yinv, Yinv, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(y, Yinv, y, prec);

    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_arb(offset, a, g);
        _arb_vec_add(offset, offset, y, g, prec);
        arb_mat_vector_mul_col(offset, cho, offset, prec);
        acb_theta_ql_sqr_dist(&dist[a], offset, cho, prec);
    }

    arb_mat_clear(Yinv);
    arb_mat_clear(cho);
    _arb_vec_clear(y, g);
    _arb_vec_clear(offset, g);
    arb_clear(pi);
}
