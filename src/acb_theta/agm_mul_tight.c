/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static void
acb_theta_agm_rel_mag_err(arf_t m, arf_t eps, acb_srcptr a, arb_srcptr d,
    slong nb, slong prec)
{
    arb_t y;
    arf_t abs;
    slong k;

    arb_init(y);
    arf_init(abs);

    arf_zero(m);
    arf_zero(eps);

    for (k = 0; k < nb; k++)
    {
        arb_zero(y);
        arb_get_ubound_arf(arb_midref(y), &d[k], prec);
        arb_exp(y, y, prec);
        arb_get_ubound_arf(abs, y, prec);
        arb_set_arf(y, abs);

        acb_get_abs_ubound_arf(abs, &a[k], prec);
        arf_mul(abs, abs, arb_midref(y), prec, ARF_RND_CEIL);
        arf_max(m, m, abs);

        acb_get_rad_ubound_arf(abs, &a[k], prec);
        arf_mul(abs, abs, arb_midref(y), prec, ARF_RND_CEIL);
        arf_max(eps, eps, abs);
    }

    arb_clear(y);
    arf_clear(abs);
}

/* This is assuming a0 corresponds to theta constants */
void
acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, int all, slong prec)
{
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong hprec = prec;
    acb_ptr v0, v;
    arf_t m0, m, eps0, eps, e;
    arb_t err;
    slong k, b;

    v0 = _acb_vec_init(n);
    v = _acb_vec_init(n);
    arf_init(m0);
    arf_init(m);
    arf_init(eps0);
    arf_init(eps);
    arf_init(e);
    arb_init(err);

    acb_theta_agm_rel_mag_err(m0, eps0, a0, d0, n, lp);
    acb_theta_agm_rel_mag_err(m, eps, a, d, n, lp);

    for (k = 0; k < n; k++)
    {
        hprec = FLINT_MAX(hprec, prec + acb_theta_sum_addprec(&d[k]));
        acb_get_mid(&v0[k], &a0[k]);
        acb_get_mid(&v[k], &a[k]);
    }

    /* Perform agm_mul or agm_sqr at high precision */
    if (a0 == a)
    {
        acb_theta_agm_mul(res, v0, v0, g, all, hprec);
    }
    else
    {
        acb_theta_agm_mul(res, v0, v, g, all, hprec);
    }

    /* New relative error wrt distances is 2^g (m0 eps + m eps0 + eps0 eps) */
    arf_mul(e, m0, eps, lp, ARF_RND_CEIL);
    arf_addmul(e, m, eps0, lp, ARF_RND_CEIL);
    arf_addmul(e, eps, eps0, lp, ARF_RND_CEIL);
    arf_mul_2exp_si(e, e, g);

    for (k = 0; k < n; k++)
    {
        arb_neg(err, &d[k]);
        arb_exp(err, err, lp);
        arb_mul_arf(err, err, e, lp);

        if (all)
        {
            for (b = 0; b < n; b++)
            {
                acb_add_error_arb(&res[k * n + b], err);
            }
        }
        else
        {
            acb_add_error_arb(&res[k], err);
        }
    }

    _acb_vec_clear(v0, n);
    _acb_vec_clear(v, n);
    arf_clear(m0);
    arf_clear(m);
    arf_clear(eps0);
    arf_clear(eps);
    arf_clear(e);
    arb_clear(err);
}
