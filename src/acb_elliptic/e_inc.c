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
/* s*(RF(x,y,1) - m*s^2*RD(x,y,1)/3), x = c^2, y=1-m*s^2 */
void
acb_elliptic_e_reduced(acb_t r, const acb_t z, const acb_t m, int times_pi, slong prec)
{
    acb_t s, c, x, y, rf, rd;

    acb_init(s);
    acb_init(c);
    acb_init(x);
    acb_init(y);
    acb_init(rf);
    acb_init(rd);

    if (times_pi)
        acb_sin_cos_pi(s, c, z, prec);
    else
        acb_sin_cos(s, c, z, prec);

    acb_mul(x, c, c, prec);
    acb_mul(y, s, s, prec);
    acb_mul(y, y, m, prec);
    acb_sub_ui(y, y, 1, prec);
    acb_neg(y, y);
    acb_one(rf);
    acb_one(rd);

    acb_elliptic_rf(rf, x, y, rf, 0, prec);
    acb_elliptic_rj(rd, x, y, rd, rd, 0, prec);

    acb_mul(y, s, s, prec);
    acb_mul(y, y, m, prec);
    acb_mul(rd, rd, y, prec);
    acb_div_ui(rd, rd, 3, prec);

    acb_sub(r, rf, rd, prec);
    acb_mul(r, r, s, prec);

    acb_clear(s);
    acb_clear(c);
    acb_clear(x);
    acb_clear(y);
    acb_clear(rf);
    acb_clear(rd);
}

void
acb_elliptic_e_inc(acb_t res, const acb_t phi, const acb_t m, int times_pi, slong prec)
{
    arb_t x, d, pi;
    acb_t z, w, r;

    if (!acb_is_finite(phi) || !acb_is_finite(m))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(m))
    {
        if (times_pi)
        {
            arb_init(pi);
            arb_const_pi(pi, prec);
            acb_mul_arb(res, phi, pi, prec);
            arb_clear(pi);
        }
        else
        {
            acb_set_round(res, phi, prec);
        }
        return;
    }

    if (acb_is_zero(phi))
    {
        acb_zero(res);
        return;
    }

    if (times_pi && acb_is_int_2exp_si(phi, -1))
    {
        acb_t t;
        acb_init(t);
        acb_mul_2exp_si(t, phi, 1);
        acb_elliptic_e(res, m, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
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

        arb_sub_ui(x, acb_realref(m), 1, prec);
        is_real = acb_is_real(phi) && acb_is_real(m) && arb_is_negative(x);

        /* left d */
        acb_zero(z);
        arf_set_mag(arb_midref(acb_realref(z)), arb_radref(d));
        mag_zero(arb_radref(d));
        arb_sub(d, d, acb_realref(z), 2 * prec + 100); /* meant to be exact */
        arb_floor(d, d, prec);

        /* w = 2 E(m) */
        acb_elliptic_e(w, m, prec);
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

        acb_elliptic_e_reduced(r, z, m, times_pi, prec);

        acb_addmul_arb(r, w, d, prec);

        /* z = phi - (d + 1) * pi */
        if (times_pi)
            acb_sub_ui(z, z, 1, prec);
        else
            acb_sub_arb(z, z, pi, prec);

        acb_elliptic_e_reduced(r2, z, m, times_pi, prec);

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

            /* w = 2 d E(m) */
            acb_elliptic_e(w, m, prec);
            acb_mul_arb(w, w, d, prec);
            acb_mul_2exp_si(w, w, 1);
        }

        acb_elliptic_e_reduced(r, z, m, times_pi, prec);
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

