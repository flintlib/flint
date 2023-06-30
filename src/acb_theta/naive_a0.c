/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_naive_a0(acb_t th, acb_srcptr z, slong nb_z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr v;
    slong k;

    v = _acb_vec_init(nb_z * n * n);
    
    acb_theta_naive_all(v, z, nb_z, tau, prec);
    for (k = 0; k < nb_z; k++)
    {
        acb_theta_get_a0(th + k * n, v + k * n * n, g);
    }
    
    _acb_vec_clear(v, nb_z * n * n);
}
