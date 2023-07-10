/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
worker_dim0(acb_ptr th, const acb_t term, slong* coords, slong g, ulong ab,
    slong ord, slong prec, slong fullprec)
{
    acb_t x;

    acb_init(x);
    acb_mul_powi(x, term, acb_theta_char_dot_slong(ab, coords, g));
    acb_add(th, th, x, fullprec);
    acb_clear(x);
}

void
acb_theta_naive_ind(acb_t th, ulong ab, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_struct* eps;
    acb_ptr c;
    acb_ptr new_z;
    int all = 0;
    slong ord = 0;
    slong nb = 1;
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, nb_z, g);
    eps = flint_malloc(nb_z * sizeof(arf_struct));
    for (k = 0; k < nb_z; k++)
    {
        arf_init(&eps[k]);
    }
    c = _acb_vec_init(nb_z);
    new_z = _acb_vec_init(g * nb_z);

    acb_theta_naive_ellipsoid(E, eps, c, new_z, ab, all, ord, z, nb_z, tau, prec);
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, new_z, tau, E, prec);

    for (k = 0; k < nb_z; k++)
    {
        acb_theta_naive_worker(&th[k], nb, &c[k], &eps[k], E, D, k, ab, ord,
            prec, worker_dim0);
    }

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    for (k = 0; k < nb_z; k++)
    {
        arf_clear(&eps[k]);        
    }
    flint_free(eps);
    _acb_vec_clear(c, nb_z);
    _acb_vec_clear(new_z, g * nb_z);
}
