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
acb_theta_dupl_transform_const_radius(arf_t rho, const arf_t r, acb_srcptr th,
    const fmpz_mat_t mat, slong prec)
{
    acb_ptr th_dupl;
    slong g = sp2gz_dim(mat);
    slong n = 1 << g;

    th_dupl = _acb_vec_init(n * n);

    acb_theta_dupl_all_const(th_dupl, th, g, prec);
    acb_theta_transform_radius(rho, r, th_dupl, mat, prec);
    acb_theta_dupl_radius(rho, rho, th, n, prec);

    _acb_vec_clear(th_dupl, n * n);
}
