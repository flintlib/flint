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
acb_theta_renormalize_const_sqr(acb_t scal, acb_srcptr th2,
                                const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = 1 + acb_theta_agm_nb_bad_steps(tau, lowprec);
    acb_ptr roots;
    acb_ptr a;
    slong n = 1 << g;

    roots = _acb_vec_init(n * nb_bad);
    a = _acb_vec_init(n);

    /* Compute lowprec roots */
    acb_theta_agm_roots(roots, tau, nb_bad, lowprec);

    /* Renormalize lowprec roots */
    acb_sqrt(scal, &th2[0], 2 * lowprec);
    acb_div(scal, scal, &roots[0], lowprec);
    _acb_vec_scalar_mul(roots, roots, n * nb_bad, scal, lowprec);

    /* Inverse agm */
    acb_theta_agm(scal, th2, roots, nb_bad, g, prec);
    acb_inv(scal, scal, prec);

    _acb_vec_clear(roots, n * nb_bad);
    _acb_vec_clear(a, n);
}
