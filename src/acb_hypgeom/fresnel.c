/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

/*
We compute the following normalized versions internally:

S(z) = (8/sqrt(pi)) int_0^z sin(2t^2) dt
C(z) = (8/sqrt(pi)) int_0^z cos(2t^2) dt

The benefit is that z^2 can be computed exactly inside erf when we have
multiplied by 1+i instead of (1+i)/sqrt(2), so we get faster evaluation
and better error bounds for Fresnel integrals on the real line (this is a
bit of a hack, and it would be better to somehow pass z^2 directly to the erf
evaluation code).
*/

void
acb_hypgeom_fresnel_erf(acb_t res1, acb_t res2, const acb_t z, slong prec)
{
    acb_t t, u, v, w1, w2;

    acb_init(t);
    acb_init(v);
    acb_init(w1);

    if (arb_is_zero(acb_imagref(z)))
    {
        acb_mul_onei(t, z);
        acb_add(w1, z, t, 2 * prec);
        acb_hypgeom_erf(t, w1, prec + 4);
        acb_mul_2exp_si(t, t, 1);

        acb_mul_onei(v, t);
        acb_add(t, t, v, prec);

        if (res1 != NULL) acb_set_arb(res1, acb_realref(t));
        if (res2 != NULL) acb_set_arb(res2, acb_imagref(t));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        acb_mul_onei(t, z);
        acb_sub(w1, t, z, 2 * prec);
        acb_hypgeom_erf(t, w1, prec + 4);
        acb_mul_2exp_si(t, t, 1);

        acb_mul_onei(v, t);
        acb_add(t, t, v, prec);

        if (res1 != NULL) acb_set_arb(res1, acb_realref(t));
        if (res1 != NULL) acb_mul_onei(res1, res1);
        if (res2 != NULL) acb_set_arb(res2, acb_imagref(t));
        if (res2 != NULL) acb_div_onei(res2, res2);
    }
    else
    {
        acb_init(u);
        acb_init(w2);

        /* w1 = (1+i)z, w2 = (1-i)z */
        acb_mul_onei(t, z);
        acb_add(w1, z, t, 2 * prec);
        acb_sub(w2, z, t, 2 * prec);

        acb_hypgeom_erf(t, w1, prec + 4);
        acb_hypgeom_erf(u, w2, prec + 4);

        /* S = (1+i) (t - ui) = (1+i) t + (1-i) u */
        /* C = (1-i) (t + ui) = (1-i) t + (1+i) u */

        acb_mul_onei(v, t);
        if (res1 != NULL) acb_add(res1, t, v, prec);
        if (res2 != NULL) acb_sub(res2, t, v, prec);

        acb_mul_onei(v, u);
        if (res1 != NULL) acb_add(res1, res1, u, prec);
        if (res1 != NULL) acb_sub(res1, res1, v, prec);
        if (res2 != NULL) acb_add(res2, res2, u, prec);
        if (res2 != NULL) acb_add(res2, res2, v, prec);

        acb_clear(u);
        acb_clear(w2);
    }

    acb_clear(t);
    acb_clear(v);
    acb_clear(w1);
}

/* derivatives: |8/sqrt(pi) sin(2z^2)|, |8/sqrt(pi) cos(2z^2)| <= 5 exp(4|xy|) */
void
acb_hypgeom_fresnel_erf_error(acb_t res1, acb_t res2, const acb_t z, slong prec)
{
    mag_t re;
    mag_t im;
    acb_t zmid;

    mag_init(re);
    mag_init(im);
    acb_init(zmid);

    if (arf_cmpabs_ui(arb_midref(acb_realref(z)), 1000) < 0 &&
        arf_cmpabs_ui(arb_midref(acb_imagref(z)), 1000) < 0)
    {
        arb_get_mag(re, acb_realref(z));
        arb_get_mag(im, acb_imagref(z));
        mag_mul(re, re, im);
        mag_mul_2exp_si(re, re, 2);
        mag_exp(re, re);
        mag_mul_ui(re, re, 5);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_mul(t, acb_realref(z), acb_imagref(z), prec);
        arb_abs(t, t);
        arb_mul_2exp_si(t, t, 2);
        arb_exp(t, t, prec);
        arb_get_mag(re, t);
        mag_mul_ui(re, re, 5);
        arb_clear(t);
    }

    mag_hypot(im, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
    mag_mul(re, re, im);

    if (arb_is_zero(acb_imagref(z)))
    {
        mag_set_ui(im, 8);  /* For real x, |S(x)| < 4, |C(x)| < 4. */
        mag_min(re, re, im);
        mag_zero(im);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        mag_set_ui(im, 8);
        mag_min(im, re, im);
        mag_zero(re);
    }
    else
    {
        mag_set(im, re);
    }

    arf_set(arb_midref(acb_realref(zmid)), arb_midref(acb_realref(z)));
    arf_set(arb_midref(acb_imagref(zmid)), arb_midref(acb_imagref(z)));

    acb_hypgeom_fresnel_erf(res1, res2, zmid, prec);

    if (res1 != NULL)
    {
        arb_add_error_mag(acb_realref(res1), re);
        arb_add_error_mag(acb_imagref(res1), im);
    }

    if (res2 != NULL)
    {
        arb_add_error_mag(acb_realref(res2), re);
        arb_add_error_mag(acb_imagref(res2), im);
    }

    mag_clear(re);
    mag_clear(im);
    acb_clear(zmid);
}

void
acb_hypgeom_fresnel(acb_t res1, acb_t res2, const acb_t z, int normalized, slong prec)
{
    slong wp;
    acb_t w;
    arb_t c;

    if (!acb_is_finite(z))
    {
        if (res1 != NULL) acb_indeterminate(res1);
        if (res2 != NULL) acb_indeterminate(res2);
        return;
    }

    acb_init(w);
    arb_init(c);

    wp = prec + 8;

    if (normalized)
    {
        arb_const_pi(c, wp);
        arb_sqrt(c, c, wp);
        arb_mul_2exp_si(c, c, -1);
        acb_mul_arb(w, z, c, wp);
        acb_hypgeom_fresnel_erf_error(res1, res2, w, wp);
    }
    else
    {
        arb_sqrt_ui(c, 2, wp);
        arb_mul_2exp_si(c, c, -1);
        acb_mul_arb(w, z, c, wp);
        acb_hypgeom_fresnel_erf_error(res1, res2, w, wp);
        arb_const_pi(c, wp);
        arb_mul_2exp_si(c, c, -1);
        arb_sqrt(c, c, wp);

        if (res1 != NULL) acb_mul_arb(res1, res1, c, wp);
        if (res2 != NULL) acb_mul_arb(res2, res2, c, wp);
    }

    if (res1 != NULL)
    {
        acb_mul_2exp_si(res1, res1, -2);
        acb_set_round(res1, res1, prec);
    }

    if (res2 != NULL)
    {
        acb_mul_2exp_si(res2, res2, -2);
        acb_set_round(res2, res2, prec);
    }

    acb_clear(w);
    arb_clear(c);
}

