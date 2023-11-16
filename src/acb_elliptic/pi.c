/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"
#include "acb_modular.h"

/* Evaluation on -pi/2 <= re(z) <= pi/2, no aliasing. */
/* s*(RF(x,y,1) + n*s^2*RJ(x,y,1,p)/3), x = c^2, y=1-m*s^2, p=1-n*s^2 */
/* complete: c = 0, s = 1 */
void
acb_elliptic_pi_reduced(acb_t r, const acb_t n,
    const acb_t z, const acb_t m, int times_pi, slong prec)
{
    acb_t s, c, x, y, p, rf, rj;

    acb_init(s);
    acb_init(c);
    acb_init(x);
    acb_init(y);
    acb_init(p);
    acb_init(rf);
    acb_init(rj);

    if (times_pi)
        acb_sin_cos_pi(s, c, z, prec);
    else
        acb_sin_cos(s, c, z, prec);

    acb_mul(x, c, c, prec);
    acb_mul(y, s, s, prec);
    acb_mul(p, y, n, prec);
    acb_mul(y, y, m, prec);
    acb_sub_ui(y, y, 1, prec);
    acb_neg(y, y);
    acb_sub_ui(p, p, 1, prec);
    acb_neg(p, p);
    acb_one(rf);
    acb_one(rj);

    acb_elliptic_rf(rf, x, y, rf, 0, prec);
    acb_elliptic_rj(rj, x, y, rj, p, 0, prec);

    acb_mul(y, s, s, prec);
    acb_mul(y, y, n, prec);
    acb_mul(rj, rj, y, prec);
    acb_div_ui(rj, rj, 3, prec);

    acb_add(r, rf, rj, prec);
    acb_mul(r, r, s, prec);

    acb_clear(s);
    acb_clear(c);
    acb_clear(x);
    acb_clear(y);
    acb_clear(p);
    acb_clear(rf);
    acb_clear(rj);
}

void
acb_elliptic_pi(acb_t r, const acb_t n, const acb_t m, slong prec)
{
    if (acb_is_zero(n))
    {
        acb_elliptic_k(r, m, prec);
    }
    else if (acb_is_zero(m))
    {
        arb_t pi;
        arb_init(pi);
        arb_const_pi(pi, prec);
        acb_sub_ui(r, n, 1, prec);
        acb_neg(r, r);
        acb_rsqrt(r, r, prec);
        acb_mul_arb(r, r, pi, prec);
        acb_mul_2exp_si(r, r, -1);
        arb_clear(pi);
    }
    else
    {
        acb_t z;
        acb_init(z);
        acb_one(z);
        acb_mul_2exp_si(z, z, -1);
        acb_elliptic_pi_reduced(r, n, z, m, 1, prec);
        acb_clear(z);
    }
}

void
acb_elliptic_pi_inc(acb_t res, const acb_t n, const acb_t phi, const acb_t m, int times_pi, slong prec)
{
    arb_t x, d, pi;
    acb_t z, w, r;

    if (!acb_is_finite(n) || !acb_is_finite(phi) || !acb_is_finite(m))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(n))
    {
        acb_elliptic_f(res, phi, m, times_pi, prec);
        return;
    }

    if (acb_is_zero(phi) || (times_pi && acb_is_int_2exp_si(phi, -1)))
    {
        acb_t t;
        acb_init(t);
        acb_mul_2exp_si(t, phi, 1);
        acb_elliptic_pi(res, n, m, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
        return;
    }

    /* Fixme: the exact argument reduction sometimes
       does not work when re(phi) = k/2 and we end up exactly on a branch cut,
       presumably due to getting the wrong sign in R_J? Should
       investigate and find a better solution. */
    if (times_pi && !acb_is_real(phi))
    {
        acb_init(z);
        arb_init(pi);
        arb_const_pi(pi, prec);
        acb_mul_arb(z, phi, pi, prec);
        acb_elliptic_pi_inc(res, n, z, m, 0, prec);
        acb_clear(z);
        arb_clear(pi);
        return;
    }

    arb_init(x);
    arb_init(d);
    arb_init(pi);
    acb_init(z);
    acb_init(w);
    acb_init(r);

    arb_set(x, acb_realref(phi));
    arb_const_pi(pi, prec);

    if (times_pi)
        arb_set(d, x);
    else
        arb_div(d, x, pi, prec);

    arb_mul_2exp_si(d, d, 1);
    arb_add_ui(d, d, 1, prec);
    arb_mul_2exp_si(d, d, -1);

    if (mag_cmp_2exp_si(arb_radref(d), -1) >= 0)
    {
        /* may span multiple periods... don't bother */
        acb_indeterminate(res);
    }
    else if (arb_contains_int(d) && !arb_is_exact(d))  /* two adjacent d */
    {
        acb_t r2, w2;
        int is_real;

        acb_init(r2);
        acb_init(w2);

        is_real = acb_is_real(phi) && acb_is_real(m) && acb_is_real(n);
        arb_sub_ui(x, acb_realref(m), 1, prec);
        is_real = is_real && arb_is_negative(x);
        arb_sub_ui(x, acb_realref(n), 1, prec);
        is_real = is_real && arb_is_negative(x);

        /* left d */
        acb_zero(z);
        arf_set_mag(arb_midref(acb_realref(z)), arb_radref(d));
        mag_zero(arb_radref(d));
        arb_sub(d, d, acb_realref(z), 2 * prec + 100); /* meant to be exact */
        arb_floor(d, d, prec);

        /* w = 2 Pi(n, m) */
        acb_elliptic_pi(w, n, m, prec);
        acb_mul_2exp_si(w, w, 1);

        /* z = phi - d * pi */
        if (times_pi)
        {
            acb_sub_arb(z, phi, d, prec);
        }
        else
        {
            arb_mul(acb_realref(z), pi, d, prec);
            arb_sub(acb_realref(z), acb_realref(phi), acb_realref(z), prec);
            arb_set(acb_imagref(z), acb_imagref(phi));
        }

        acb_elliptic_pi_reduced(r, n, z, m, times_pi, prec);

        acb_addmul_arb(r, w, d, prec);

        /* z = phi - (d + 1) * pi */
        if (times_pi)
            acb_sub_ui(z, z, 1, prec);
        else
            acb_sub_arb(z, z, pi, prec);

        acb_elliptic_pi_reduced(r2, n, z, m, times_pi, prec);

        arb_add_ui(d, d, 1, prec);
        acb_addmul_arb(r2, w, d, prec);

        arb_union(acb_realref(res), acb_realref(r), acb_realref(r2), prec);
        arb_union(acb_imagref(res), acb_imagref(r), acb_imagref(r2), prec);

        if (is_real)
            arb_zero(acb_imagref(res));

        acb_clear(r2);
        acb_clear(w2);
    }
    else
    {
        /* this could still be inexact if d is large (which is fine) */
        arb_floor(d, d, prec);

        if (arb_is_zero(d))
        {
            acb_set(z, phi);
            acb_zero(w);
        }
        else
        {
            /* z = phi - d*pi */
            if (times_pi)
            {
                acb_sub_arb(z, phi, d, prec);
            }
            else
            {
                arb_mul(acb_realref(z), pi, d, prec);
                arb_sub(acb_realref(z), acb_realref(phi), acb_realref(z), prec);
                arb_set(acb_imagref(z), acb_imagref(phi));
            }

            /* w = 2 d Pi(n, m) */
            acb_elliptic_pi(w, n, m, prec);
            acb_mul_arb(w, w, d, prec);
            acb_mul_2exp_si(w, w, 1);
        }

        acb_elliptic_pi_reduced(r, n, z, m, times_pi, prec);

        acb_add(r, r, w, prec);
        acb_set(res, r);
    }

    arb_clear(x);
    arb_clear(d);
    arb_clear(pi);
    acb_clear(z);
    acb_clear(w);
    acb_clear(r);
}

