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
acb_theta_newton_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau,
                         slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_t scal1, scal2;

    acb_init(scal1);
    acb_init(scal2);

    acb_theta_newton_half_proj(th2, z, tau, prec);
    acb_theta_dupl_all(th2, th2, g, prec);
    acb_theta_renormalize_sqr(scal1, scal2, th2, th2 + n * n, z, tau, prec);

    _acb_vec_scalar_mul(th2, th2, n * n, scal1, prec);
    _acb_vec_scalar_mul(th2 + n * n, th2 + n * n, n * n, scal2, prec);

    acb_clear(scal1);
    acb_clear(scal2);
}
