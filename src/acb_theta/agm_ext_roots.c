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
acb_theta_agm_ext_roots(acb_ptr roots, acb_srcptr z, const acb_mat_t tau,
                        slong nb_bad, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr z_0;
    slong k;

    acb_mat_init(w, g, g);
    z_0 = _acb_vec_init(2 * g);

    acb_mat_set(w, tau);
    _acb_vec_set(z_0, z, g);

    for (k = 0; k < nb_bad; k++)
    {
        acb_theta_naive(&roots[2 * k * n], z_0, 2, w, prec);
        acb_mat_scalar_mul_2exp_si(w, w, 1);
    }

    acb_mat_clear(w);
    _acb_vec_clear(z_0, 2 * g);
}
