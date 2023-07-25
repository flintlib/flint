/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* This is assuming a0 corresponds to theta constants */
void
acb_theta_agm_mul_tight(acb_ptr r, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong hprec = prec;
    acb_ptr v0, v;
    arf_t m0, m, eps0, eps, e, t;
    arb_t err;
    slong k;

    v0 = _acb_vec_init(n);
    v = _acb_vec_init(n);
    arf_init(m0);
    arf_init(m);
    arf_init(eps0);
    arf_init(eps);
    arf_init(e);
    arf_init(t);
    arb_init(err);

    acb_theta_agm_rel_mag_err(m0, eps0, a0, d0, n, prec);
    acb_theta_agm_rel_mag_err(m, eps, a, d, n, prec);

    for (k = 0; k < n; k++)
    {
        hprec = FLINT_MAX(hprec, prec + acb_theta_dist_addprec(&d[k]));
        acb_get_mid(&v0[k], &a0[k]);
        acb_get_mid(&v[k], &a[k]);
    }

    /* Perform agm_mul or agm_sqr at high precision */
    if (a0 == a)
    {
        acb_theta_agm_sqr(r, v0, g, hprec);
    }
    else
    {
        acb_theta_agm_mul(r, v0, v, g, hprec);
    }

    /* New relative error wrt distances is m0 eps + m eps0 + eps0 eps */
    arf_mul(e, m0, eps, lp, ARF_RND_CEIL);
    arf_mul(t, m, eps0, lp, ARF_RND_CEIL);
    arf_add(e, e, t, lp, ARF_RND_CEIL);
    arf_mul(t, eps, eps0, lp, ARF_RND_CEIL);
    arf_add(e, e, t, lp, ARF_RND_CEIL);

    for (k = 0; k < n; k++)
    {
        arb_neg(err, &d[k]);
        arb_exp(err, err, prec);
        arb_mul_arf(err, err, e, lp);
        acb_add_error_arb(&r[k], err);
    }

    _acb_vec_clear(v0, n);
    _acb_vec_clear(v, n);
    arf_clear(m0);
    arf_clear(m);
    arf_clear(eps0);
    arf_clear(eps);
    arf_clear(e);
    arf_clear(t);
    arb_clear(err);
}
