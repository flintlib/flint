/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define log __builtin_log
#else
# include <math.h>
#endif

/* IMAG: erf(z) = 2z/sqrt(pi) * 1F1(1/2, 3/2, -z^2) */
void
acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, slong prec)
{
    acb_t a, t, w;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(t);
    acb_init(w);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_set_ui(b, 3);
    acb_mul_2exp_si(b, b, -1);
    acb_one(b + 1);

    acb_mul(w, z, z, prec);
    acb_neg(w, w);

    acb_hypgeom_pfq_direct(t, a, 1, b, 2, w, -1, prec);

    acb_mul(t, t, z, prec);
    arb_const_sqrt_pi(acb_realref(w), prec);
    acb_div_arb(t, t, acb_realref(w), prec);

    acb_mul_2exp_si(res, t, 1);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(t);
    acb_clear(w);
}

/* REAL: erf(x) = 2x/sqrt(pi) * exp(-x^2) 1F1(1, 3/2, x^2) */
void
acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, slong prec)
{
    acb_t a, b, t, w;

    acb_init(a);
    acb_init(b);
    acb_init(t);
    acb_init(w);

    acb_set_ui(b, 3);
    acb_mul_2exp_si(b, b, -1);

    acb_mul(w, z, z, prec);

    acb_hypgeom_pfq_direct(t, a, 0, b, 1, w, -1, prec);

    acb_neg(w, w);
    acb_exp(w, w, prec);
    acb_mul(t, t, w, prec);

    acb_mul(t, t, z, prec);
    arb_const_sqrt_pi(acb_realref(w), prec);
    acb_div_arb(t, t, acb_realref(w), prec);

    acb_mul_2exp_si(res, t, 1);

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
    acb_clear(w);
}

void
acb_hypgeom_erf_asymp(acb_t res, const acb_t z, int complementary, slong prec, slong prec2)
{
    acb_t a, t, u;

    acb_init(a);
    acb_init(t);
    acb_init(u);

    if (!acb_is_exact(z) &&
        (arf_cmpabs_ui(arb_midref(acb_realref(z)), prec) < 0) &&
        (arf_cmpabs_ui(arb_midref(acb_imagref(z)), prec) < 0))
    {
        acb_t zmid;
        mag_t re_err, im_err;

        acb_init(zmid);
        mag_init(re_err);
        mag_init(im_err);

        acb_hypgeom_erf_propagated_error(re_err, im_err, z);
        arf_set(arb_midref(acb_realref(zmid)), arb_midref(acb_realref(z)));
        arf_set(arb_midref(acb_imagref(zmid)), arb_midref(acb_imagref(z)));

        acb_hypgeom_erf_asymp(res, zmid, complementary, prec, prec2);

        arb_add_error_mag(acb_realref(res), re_err);
        arb_add_error_mag(acb_imagref(res), im_err);

        acb_clear(zmid);
        mag_clear(re_err);
        mag_clear(im_err);

        return;
    }

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_mul(t, z, z, prec2);

    acb_hypgeom_u_asymp(u, a, a, t, -1, prec2);

    acb_neg(t, t);
    acb_exp(t, t, prec2);
    acb_mul(u, u, t, prec2);

    arb_const_sqrt_pi(acb_realref(t), prec2);
    arb_zero(acb_imagref(t));
    acb_mul(t, t, z, prec2);
    acb_div(u, u, t, prec2);

    /* branch cut term: -1 or 1 */
    acb_csgn(acb_realref(t), z);
    arb_zero(acb_imagref(t));

    if (complementary)
    {
        /* erfc(z) = 1 - erf(z) = u - (sgn - 1) */
        acb_sub_ui(t, t, 1, prec);
        acb_sub(t, u, t, prec);
    }
    else
    {
        /* erf(z) = sgn - u */
        acb_sub(t, t, u, prec);
    }

    if (arb_is_zero(acb_imagref(z)))
    {
        arb_zero(acb_imagref(t));
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        if (complementary)
            arb_one(acb_realref(t));
        else
            arb_zero(acb_realref(t));
    }

    acb_set(res, t);

    acb_clear(a);
    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_erf_propagated_error(mag_t re, mag_t im, const acb_t z)
{
    mag_t x, y;

    mag_init(x);
    mag_init(y);

    /* |exp(-(x+y)^2)| = exp(y^2-x^2) */
    arb_get_mag(y, acb_imagref(z));
    mag_mul(y, y, y);

    arb_get_mag_lower(x, acb_realref(z));
    mag_mul_lower(x, x, x);

    if (mag_cmp(y, x) >= 0)
    {
        mag_sub(re, y, x);
        mag_exp(re, re);
    }
    else
    {
        mag_sub_lower(re, x, y);
        mag_expinv(re, re);
    }

    /* Radius. */
    mag_hypot(x, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
    mag_mul(re, re, x);

    /* 2/sqrt(pi) < 289/256 */
    mag_mul_ui(re, re, 289);
    mag_mul_2exp_si(re, re, -8);

    if (arb_is_zero(acb_imagref(z)))
    {
        /* todo: could bound magnitude even for complex numbers */
        mag_set_ui(y, 2);
        mag_min(re, re, y);

        mag_zero(im);
    }
    else if (arb_is_zero(acb_realref(z)))
    {
        mag_swap(im, re);
        mag_zero(re);
    }
    else
    {
        mag_set(im, re);
    }

    mag_clear(x);
    mag_clear(y);
}

void
acb_hypgeom_erf_1f1(acb_t res, const acb_t z, slong prec,
    slong wp, int more_imaginary)
{
    if (acb_rel_accuracy_bits(z) >= wp)
    {
        if (more_imaginary)
            acb_hypgeom_erf_1f1a(res, z, wp);
        else
            acb_hypgeom_erf_1f1b(res, z, wp);
    }
    else
    {
        acb_t zmid;
        mag_t re_err, im_err;

        acb_init(zmid);
        mag_init(re_err);
        mag_init(im_err);

        acb_hypgeom_erf_propagated_error(re_err, im_err, z);
        arf_set(arb_midref(acb_realref(zmid)), arb_midref(acb_realref(z)));
        arf_set(arb_midref(acb_imagref(zmid)), arb_midref(acb_imagref(z)));

        if (more_imaginary)
            acb_hypgeom_erf_1f1a(res, zmid, wp);
        else
            acb_hypgeom_erf_1f1b(res, zmid, wp);

        arb_add_error_mag(acb_realref(res), re_err);
        arb_add_error_mag(acb_imagref(res), im_err);

        acb_clear(zmid);
        mag_clear(re_err);
        mag_clear(im_err);
    }

    acb_set_round(res, res, prec);
}

void
acb_hypgeom_erf(acb_t res, const acb_t z, slong prec)
{
    double x, y, abs_z2, log_z, log_erf_z_asymp;
    slong prec2, wp;
    int more_imaginary;

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_real(z))
    {
        arb_hypgeom_erf(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (acb_is_zero(z))
    {
        acb_zero(res);
        return;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -64) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -64) < 0))
    {
        acb_hypgeom_erf_1f1(res, z, prec, prec, 1);
        return;
    }

    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        acb_hypgeom_erf_asymp(res, z, 0, prec, prec);
        return;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    abs_z2 = x * x + y * y;
    log_z = 0.5 * log(abs_z2);
    /* estimate of log(erf(z)), disregarding csgn term */
    log_erf_z_asymp = y*y - x*x - log_z;

    if (log_z - abs_z2 < -(prec + 8) * 0.69314718055994530942)
    {
        /* If the asymptotic term is small, we can
           compute with reduced precision. */
        prec2 = FLINT_MIN(prec + 4 + log_erf_z_asymp * 1.4426950408889634074, (double) prec);
        prec2 = FLINT_MAX(8, prec2);
        prec2 = FLINT_MIN(prec2, prec);

        acb_hypgeom_erf_asymp(res, z, 0, prec, prec2);
    }
    else
    {
        more_imaginary = arf_cmpabs(arb_midref(acb_imagref(z)),
                                    arb_midref(acb_realref(z))) > 0;

        /* Worst case: exp(|x|^2), computed: exp(x^2).
           (x^2+y^2) - (x^2-y^2) = 2y^2, etc. */
        if (more_imaginary)
            wp = prec + FLINT_MAX(2 * x * x, 0.0) * 1.4426950408889634074 + 5;
        else
            wp = prec + FLINT_MAX(2 * y * y, 0.0) * 1.4426950408889634074 + 5;

        acb_hypgeom_erf_1f1(res, z, prec, wp, more_imaginary);
    }
}

