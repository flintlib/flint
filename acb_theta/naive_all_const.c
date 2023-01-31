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
acb_theta_naive_all_const(acb_ptr th, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_ptr z;

    z = _acb_vec_init(g);

    _acb_vec_zero(z, g);
    acb_theta_naive_all(th, z, 1, tau, prec);

    _acb_vec_clear(z, g);
}
