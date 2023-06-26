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
acb_theta_renormalize_sqr(acb_t scal_z, acb_t scal_0, acb_srcptr th2_z,
                          acb_srcptr th2_0, acb_srcptr z, const acb_mat_t tau,
                          slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = acb_theta_agm_ext_nb_bad_steps(z, tau, lowprec);
    acb_ptr a;
    acb_ptr th2;
    acb_ptr roots;
    acb_t scal;
    slong n = 1 << g;
    slong k;

    a = _acb_vec_init(2 * n);
    th2 = _acb_vec_init(2 * n);
    roots = _acb_vec_init(2 * n * nb_bad);
    acb_init(scal);

    _acb_vec_set(th2, th2_z, n);
    _acb_vec_set(th2 + n, th2_0, n);

    /* Compute lowprec square roots */
    acb_theta_agm_ext_roots(roots, z, tau, nb_bad, prec);

    /* Renormalize lowprec square roots */
    acb_sqrt(scal, &th2[n], 2 * lowprec);
    acb_div(scal, scal, &roots[n], lowprec);
    _acb_vec_scalar_mul(roots, roots, 2 * n * nb_bad, scal, lowprec);

    acb_sqrt(scal, &th2[0], 2 * lowprec);
    acb_div(scal, scal, &roots[0], lowprec);
    for (k = 0; k < nb_bad; k++)
    {
        _acb_vec_scalar_mul(&roots[2 * n * k], &roots[2 * n * k], n, scal,
                            lowprec);
        acb_sqrt(scal, scal, lowprec);
    }

    /* Inverse agm */
    acb_theta_agm_ext(scal_z, scal_0, th2, roots, nb_bad, g, prec);
    acb_mul(scal_z, scal_z, scal_0, prec);
    acb_inv(scal_z, scal_z, prec);
    acb_inv(scal_0, scal_0, prec);

    _acb_vec_clear(a, 2 * n);
    _acb_vec_clear(th2, 2 * n);
    _acb_vec_clear(roots, 2 * n * nb_bad);
    acb_clear(scal);
}
