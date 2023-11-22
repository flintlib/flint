/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static void
acb_theta_agm_rel_mag_err(arf_t m, arf_t eps, acb_srcptr a, arb_srcptr d,
    slong nb, slong prec)
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

    for (k = 0; k < nb; k++)
    {
        arb_zero(y);
        arb_get_ubound_arf(arb_midref(y), &d[k], prec);
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

/* This is assuming a0 corresponds to theta constants */
static void
acb_theta_agm_mul_tight_gen(acb_ptr res, acb_srcptr a0, acb_srcptr a,
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
        acb_theta_agm_mul(res, v0, v0, g, hprec);
    }
    else
    {
        acb_theta_agm_mul(res, v0, v, g, hprec);
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
        acb_add_error_arb(&res[k], err);
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

/* there might be a way of saving some multiplications here by writing in terms
   of real & imaginary parts */
static void
acb_theta_agm_mul_tight_g1(acb_ptr res, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    acb_t t;
    acb_ptr aux;

    acb_init(t);
    aux = _acb_vec_init(2);

    if (a == a0)
    {
        acb_sqr(t, &a[0], prec);
        acb_sqr(&aux[0], &a[1], prec);
        acb_add(&aux[0], &aux[0], t, prec + acb_theta_dist_addprec(&d[0]));
        acb_mul(&aux[1], &a[0], &a[1], prec);
        acb_mul_2exp_si(&aux[1], &aux[1], 1);
    }
    else
    {
        acb_mul(t, &a0[0], &a[0], prec);
        acb_mul(&aux[0], &a0[1], &a[1], prec);
        acb_add(&aux[0], &aux[0], t, prec + acb_theta_dist_addprec(&d[0]));
        acb_mul(t, &a0[0], &a[1], prec);
        acb_mul(&aux[1], &a0[1], &a[0], prec);
        acb_add(&aux[1], &aux[1], t, prec + acb_theta_dist_addprec(&d[1]));
    }
    _acb_vec_scalar_mul_2exp_si(res, aux, 2, -1);

    acb_clear(t);
    _acb_vec_clear(aux, 2);
}

void
acb_theta_agm_mul_tight(acb_ptr res, acb_srcptr a0, acb_srcptr a,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    if (g == 1)
    {
        acb_theta_agm_mul_tight_g1(res, a0, a, d0, d, g, prec);
    }
    else
    {
        acb_theta_agm_mul_tight_gen(res, a0, a, d0, d, g, prec);
    }
}
