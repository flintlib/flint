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
acb_theta_newton_const_sqr(acb_ptr th2, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_t scal;

    acb_init(scal);

    acb_theta_newton_const_half_proj(th2, tau, prec);
    acb_theta_dupl_const(th2, th2, g, prec);
    acb_theta_renormalize_const_sqr(scal, th2, tau, prec);
    _acb_vec_scalar_mul(th2, th2, 1 << g, scal, prec);

    acb_clear(scal);
}
