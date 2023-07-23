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
acb_theta_agm_rel_mag_err(arf_t m, arf_t eps, acb_srcptr a, arb_srcptr dist,
    slong n, slong prec)
{
    acb_t x, err;
    arb_t y;
    arf_t abs;
    slong k;

    acb_init(x);
    acb_init(err);
    arb_init(y);
    arf_init(abs);

    arf_zero(m);
    arf_zero(eps);
    
    for (k = 0; k < n; k++)
    {
        arb_neg(y, &dist[k]);
        arb_exp(y, y, prec);
        acb_mul_arb(x, &a[k], y, prec);

        acb_abs(y, x, prec);
        arb_get_ubound_arf(abs, y, prec);
        arf_max(m, m, abs);

        acb_zero(err);
        arf_set_mag(arb_midref(acb_realref(err)), arb_radref(acb_realref(x)));
        arf_set_mag(arb_midref(acb_imagref(err)), arb_radref(acb_imagref(x)));
        acb_abs(y, err, prec);
        arb_get_ubound_arf(abs, y, prec);
        arf_max(eps, eps, abs);
    }

    acb_clear(x);
    acb_clear(err);
    arb_clear(y);
    arf_clear(abs);
}

/* This is assuming a1 corresponds to theta constants */
void
acb_theta_agm_mul_tight(acb_ptr r, acb_srcptr a1, acb_srcptr a2,
    arb_srcptr d1, arb_srcptr d2, slong g, slong prec)
{
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong hprec = prec;
    acb_ptr v1, v2;
    arf_t m1, m2, eps1, eps2, eps, t;
    arb_t err;
    ulong a;

    v1 = _acb_vec_init(n);
    v2 = _acb_vec_init(n);
    arf_init(m1);
    arf_init(m2);
    arf_init(eps1);
    arf_init(eps2);
    arf_init(eps);
    arf_init(t);
    arb_init(err);

    acb_theta_agm_rel_mag_err(m1, eps1, a1, d1, n, lp);
    acb_theta_agm_rel_mag_err(m2, eps2, a2, d2, n, lp);
    for (a = 0; a < n; a++)
    {
        hprec = FLINT_MAX(hprec, prec + acb_theta_dist_addprec(&d2[a]));
        acb_get_mid(&v1[a], &a1[a]);
        acb_get_mid(&v2[a], &a2[a]);
    }

    /* Perform agm_mul or agm_sqr at high precision */
    if (a1 == a2)
    {
        acb_theta_agm_sqr(r, v1, g, hprec);
    }
    else
    {
        acb_theta_agm_mul(r, v1, v2, g, hprec);
    }

    /* New relative error wrt distances is m1 eps2 + m2 eps1 + eps1 eps2 */
    arf_mul(eps, m1, eps2, lp, ARF_RND_CEIL);
    arf_mul(t, m2, eps1, lp, ARF_RND_CEIL);
    arf_add(eps, eps, t, lp, ARF_RND_CEIL);
    arf_mul(t, eps2, eps1, lp, ARF_RND_CEIL);
    arf_add(eps, eps, t, lp, ARF_RND_CEIL);
    
    for (a = 0; a < n; a++)
    {
        arb_mul_arf(err, &d2[a], eps, lp);
        acb_add_error_arb(&r[a], err);
    }

    _acb_vec_clear(v1, g);
    _acb_vec_clear(v2, g);
    arf_clear(m1);
    arf_clear(m2);
    arf_clear(eps1);
    arf_clear(eps2);
    arf_clear(eps);
    arf_clear(t);
    arb_clear(err);    
}
