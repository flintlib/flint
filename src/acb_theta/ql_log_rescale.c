/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_ql_log_rescale(acb_t f, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t Yinv;
    arb_ptr y;

    arb_mat_init(Yinv, g, g);
    y = _arb_vec_init(g);

    acb_mat_get_imag(Yinv, tau);
    arb_mat_inv(Yinv, Yinv, prec);
    _acb_vec_get_imag(y, z, g);

    acb_zero(f);
    arb_mat_bilinear_form(acb_imagref(f), Yinv, y, y, prec);

    arb_mat_clear(Yinv);
    _arb_vec_clear(y, g);
}
