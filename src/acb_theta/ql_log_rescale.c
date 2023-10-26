/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_ql_log_rescale(acb_t res, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t Yinv;
    arb_ptr y, w;
    int b;

    arb_mat_init(Yinv, g, g);
    y = _arb_vec_init(g);
    w = _arb_vec_init(g);

    acb_mat_get_imag(Yinv, tau);
    b = arb_mat_inv(Yinv, Yinv, prec);
    if (!b)
    {
        arb_mat_indeterminate(Yinv);
    }
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(w, Yinv, y, prec);

    acb_zero(res);
    arb_dot(acb_imagref(res), NULL, 0, y, 1, w, 1, g, prec);

    arb_mat_clear(Yinv);
    _arb_vec_clear(y, g);
    _arb_vec_clear(w, g);
}
