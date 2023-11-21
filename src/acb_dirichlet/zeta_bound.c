/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_mag_pow(mag_t res, const mag_t x, const mag_t y)
{
    arb_t t, u;
    arb_init(t);
    arb_init(u);

    arf_set_mag(arb_midref(t), x);
    arf_set_mag(arb_midref(u), y);
    arb_pow(t, t, u, 2 * MAG_BITS);
    arb_get_mag(res, t);

    arb_clear(t);
    arb_clear(u);
}

/* zeta(1+s) < 1+1/s */
static void
mag_zeta1p(mag_t res, const mag_t s)
{
    mag_t t;
    mag_init(t);
    mag_one(t);
    mag_div(t, t, s);
    mag_add_ui(res, t, 1);
    mag_clear(t);
}

/* |zeta(s)| <= (2pi)^sigma |gamma(1-s)| exp(pi|t|/2) zeta(1-sigma) / pi */
void
acb_dirichlet_zeta_bound_functional_equation(mag_t res, const acb_t s)
{
    slong prec, p;
    acb_t z;
    arb_t x;
    mag_t t;

    if (!arb_is_negative(acb_realref(s)))
    {
        mag_inf(res);
        return;
    }

    acb_init(z);
    arb_init(x);
    mag_init(t);

    prec = 0;
    p = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(s)));
    prec = FLINT_MAX(p, prec);
    p = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(s)));
    prec = FLINT_MAX(p, prec);
    /* we *could* increase the precision even further to return finite
       results for huge input... but there's not much practical reason to... */
    prec = FLINT_MIN(prec, 1000);

    prec += MAG_BITS;

    /* gamma(1-s) */
    acb_sub_ui(z, s, 1, prec);
    acb_neg(z, z);
    acb_gamma(z, z, prec);
    acb_get_mag(res, z);

    /* (2pi)^sigma */
    arb_const_pi(x, prec);
    arb_mul_2exp_si(x, x, 1);
    arb_pow(x, x, acb_realref(s), prec);
    arb_get_mag(t, x);
    mag_mul(res, res, t);

    /* 1/pi */
    mag_div_ui(res, res, 3);

    /* exp(pi|t|/2) */
    arb_const_pi(x, prec);
    arb_mul(x, x, acb_imagref(s), prec);
    arb_abs(x, x);
    arb_mul_2exp_si(x, x, -1);
    arb_exp(x, x, prec);
    arb_get_mag(t, x);
    mag_mul(res, res, t);

    /* zeta(1-s) */
    arb_neg(x, acb_realref(s));
    arb_get_mag_lower(t, x);
    mag_zeta1p(t, t);
    mag_mul(res, res, t);

    acb_clear(z);
    arb_clear(x);
    mag_clear(t);
}

/*
Rademacher 43.3:
Assume -eta <= sigma <= 1 + eta where 0 < eta <= 1/2. Then:

|zeta(s)| < 3 |(1+s)/(1-s)| |(1+s)/(2pi)|^e zeta(1+eta)

where e = (1+eta-sigma)/2. Inside the strip, we use this formula with
eta = 0.1 (this could be improved).
*/
void
acb_dirichlet_zeta_bound_strip(mag_t res, const acb_t s)
{
    arf_t eta, a;
    acb_t s1;
    mag_t t, u, v;

    acb_init(s1);
    arf_init(eta);
    arf_init(a);
    mag_init(t);
    mag_init(u);
    mag_init(v);

    /* We need -eta <= sigma <= 1 + eta where sigma = m +/- r,
       i.e. eta >= max(-m, m-1) + r. */
    arf_neg(eta, arb_midref(acb_realref(s)));
    arf_sub_ui(a, arb_midref(acb_realref(s)), 1, MAG_BITS, ARF_RND_CEIL);
    arf_max(eta, eta, a);
    arf_set_mag(a, arb_radref(acb_realref(s)));
    arf_add(eta, eta, a, MAG_BITS, ARF_RND_CEIL);
    /* eta = max(eta, 0.1), to avoid the pole */
    arf_set_d(a, 0.1);
    arf_max(eta, eta, a);

    /* Requires 0 <= eta <= 1/2. */
    if (arf_cmpabs_2exp_si(eta, -1) <= 0)
    {
        /* t = |1+s|/(2pi) */
        acb_add_ui(s1, s, 1, MAG_BITS);
        acb_get_mag(t, s1);
        mag_set_ui_2exp_si(u, 163, -10); /* 1/(2pi) < 163/1024 */
        mag_mul(t, t, u);

        /* a = (1+eta-sigma)/2 */
        arf_set_mag(a, arb_radref(acb_realref(s)));
        arf_add(a, eta, a, MAG_BITS, ARF_RND_CEIL);
        arf_sub(a, a, arb_midref(acb_realref(s)), MAG_BITS, ARF_RND_CEIL);
        arf_add_ui(a, a, 1, MAG_BITS, ARF_RND_CEIL);
        arf_mul_2exp_si(a, a, -1);
        if (arf_sgn(a) < 0)
            arf_zero(a);
        arf_get_mag(u, a);

        /* t = (|1+s|/(2pi))^((1+eta-sigma)/2) */
        _mag_pow(t, t, u);

        /* 3|1+s|/|1-s| */
        acb_get_mag(u, s1);
        mag_mul(t, t, u);
        acb_sub_ui(s1, s, 1, MAG_BITS);
        acb_get_mag_lower(u, s1);
        mag_div(t, t, u);
        mag_mul_ui(t, t, 3);

        /* zeta(1+eta) */
        arf_get_mag_lower(u, eta);
        mag_zeta1p(u, u);
        mag_mul(t, t, u);

        mag_set(res, t);
    }
    else
    {
        mag_inf(res);
    }

    acb_clear(s1);
    arf_clear(eta);
    arf_clear(a);
    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
}

/*
We have three cases: Rademacher's bound valid on -0.5 <= sigma <= 1.5, the
trivial bound on sigma > 1, and the functional equation valid on sigma < 0.
We intersect s with three separate domains for evaluation: right, inside,
and left of the extended strip [-0.25,1.25].

Sharper bounds could be used precisely when sigma = 1/2,
or very close to the line sigma = 1.
*/
void
acb_dirichlet_zeta_bound(mag_t res, const acb_t s)
{
    arb_t strip;
    mag_t t;

    if (!acb_is_finite(s))
    {
        mag_inf(res);
        return;
    }

    arb_init(strip);
    mag_init(t);

    arf_set_ui_2exp_si(arb_midref(strip), 1, -1);
    mag_set_ui_2exp_si(arb_radref(strip), 3, -2);

    if (arb_le(strip, acb_realref(s)))
    {
        arb_get_mag_lower(res, acb_realref(s));
        mag_one(t);
        mag_sub_lower(res, res, t);
        mag_zeta1p(res, res);
    }
    else if (arb_contains(strip, acb_realref(s)))
    {
        acb_dirichlet_zeta_bound_strip(res, s);
    }
    else if (arb_le(acb_realref(s), strip))
    {
        acb_dirichlet_zeta_bound_functional_equation(res, s);
    }
    else
    {
        acb_t ss;
        arf_t x1, x2;

        acb_init(ss);
        arf_init(x1);
        arf_init(x2);

        /* Since s overlaps at least two regions, it must certainly
           overlap with the extended strip. */
        arb_set(acb_imagref(ss), acb_imagref(s));
        arb_intersection(acb_realref(ss), acb_realref(s), strip, MAG_BITS);
        acb_dirichlet_zeta_bound_strip(res, ss);

        /* We may have real parts > 1.25. */
        /* The bound computed for the extended strip *should* already be
           larger than zeta(1.25) < 5, but just to be sure... */
        mag_set_ui(t, 5);
        mag_max(res, res, t);

        /* Finally, we may have have real parts < -0.25. */
        arf_set_mag(x1, arb_radref(acb_realref(s)));
        arf_sub(x1, arb_midref(acb_realref(s)), x1, MAG_BITS, ARF_RND_FLOOR);
        arf_set_d(x2, -0.25);
        if (arf_cmp(x1, x2) < 0)
        {
            arb_set_interval_arf(acb_realref(ss), x1, x2, MAG_BITS);
            acb_dirichlet_zeta_bound_functional_equation(t, ss);
            mag_max(res, res, t);
        }

        acb_clear(ss);
        arf_clear(x1);
        arf_clear(x2);
    }

    arb_clear(strip);
    mag_clear(t);
}

/*
  |f'(s)|  <= |f(s +/- R)| / R
  |f''(s)| <= 2 |f(s +/- R)| / R^2
*/
void
acb_dirichlet_zeta_deriv_bound(mag_t der1, mag_t der2, const acb_t s)
{
    mag_t R, M;
    acb_t t;

    mag_init(R);
    mag_init(M);
    acb_init(t);

    /* R = 1/8 */
    mag_set_ui_2exp_si(R, 1, -3);

    /* t = s +/- R */
    acb_set(t, s);
    mag_add(arb_radref(acb_realref(t)), arb_radref(acb_realref(t)), R);
    mag_add(arb_radref(acb_imagref(t)), arb_radref(acb_imagref(t)), R);
    /* M = |f(s +/- R)| */
    acb_dirichlet_zeta_bound(M, t);
    /* der1 = |f'(s)| */
    mag_div(der1, M, R);
    /* der2 = |f''(s)| */
    mag_div(der2, der1, R);
    mag_mul_2exp_si(der2, der2, 1);

    acb_clear(t);
    mag_clear(R);
    mag_clear(M);
}

