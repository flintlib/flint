/*
    Copyright (C) 2015, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_si_asymp(acb_t res, const acb_t z, slong prec)
{
    acb_t t, u, w, v, one;

    acb_init(t);
    acb_init(u);
    acb_init(w);
    acb_init(v);
    acb_init(one);

    acb_one(one);
    acb_mul_onei(w, z);

    /* u = U(1,1,iz) */
    acb_hypgeom_u_asymp(u, one, one, w, -1, prec);
    /* v = e^(-iz) */
    acb_neg(v, w);
    acb_exp(v, v, prec);
    acb_mul(t, u, v, prec);

    if (acb_is_real(z))
    {
        arb_div(acb_realref(t), acb_realref(t), acb_realref(z), prec);
        arb_zero(acb_imagref(t));
        acb_neg(t, t);
    }
    else
    {
        /* u = U(1,1,-iz) */
        acb_neg(w, w);
        acb_hypgeom_u_asymp(u, one, one, w, -1, prec);
        acb_inv(v, v, prec);
        acb_addmul(t, u, v, prec);

        acb_div(t, t, z, prec);
        acb_mul_2exp_si(t, t, -1);
        acb_neg(t, t);
    }

    if (arb_is_zero(acb_realref(z)))
    {
        /* the function is imaginary */
        arb_zero(acb_realref(t));
    }
    else if (arb_is_positive(acb_realref(z)))
    {
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, -1);
        arb_add(acb_realref(t), acb_realref(t), acb_realref(u), prec);
    }
    else if (arb_is_negative(acb_realref(z)))
    {
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, -1);
        arb_sub(acb_realref(t), acb_realref(t), acb_realref(u), prec);
    }
    else
    {
        /* add [-pi,pi]/2 */
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, -1);
        arb_add_error(acb_imagref(t), acb_realref(u));
    }

    acb_swap(res, t);

    if (!acb_is_finite(res))
        acb_indeterminate(res);

    acb_clear(t);
    acb_clear(u);
    acb_clear(w);
    acb_clear(v);
    acb_clear(one);
}

void
acb_hypgeom_si_1f2(acb_t res, const acb_t z, slong prec)
{
    acb_t a, t;
    acb_struct b[3];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(b + 2);
    acb_init(t);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_set_ui(b, 3);
    acb_mul_2exp_si(b, b, -1);
    acb_set(b + 1, b);
    acb_one(b + 2);

    acb_mul(t, z, z, prec);
    acb_mul_2exp_si(t, t, -2);
    acb_neg(t, t);
    acb_hypgeom_pfq_direct(t, a, 1, b, 3, t, -1, prec);
    acb_mul(t, t, z, prec);

    acb_swap(res, t);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(b + 2);
    acb_clear(t);
}

/* Bound propagated error |si(z)-si(mid(z))| using
   |si'(z)| = |sin(z)/z| <= cosh(|im(z)|)/max(1, |z|).
*/
static void
acb_hypgeom_si_prop_error(mag_t err, const acb_t z)
{
    if (acb_is_exact(z))
    {
        mag_zero(err);
    }
    else
    {
        mag_t t, u;

        mag_init(t);
        mag_init(u);

        arb_get_mag(t, acb_imagref(z));
        mag_cosh(t, t);
        acb_get_mag_lower(u, z);
        if (mag_cmp_2exp_si(u, 0) < 0)
            mag_one(u);
        mag_div(t, t, u);
        mag_hypot(u, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul(err, t, u);

        mag_clear(t);
        mag_clear(u);
    }
}

void
acb_hypgeom_si(acb_t res, const acb_t z, slong prec)
{
    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
    }
    else if (acb_is_real(z))
    {
        arb_hypgeom_si(acb_realref(res), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        double asymp_accuracy, a, b, absz, cancellation, rlog2 = 1.4426950408889634;
        slong wp;
        int use_asymp;
        acb_t m;
        mag_t err;
        int pure_imag = arb_is_zero(acb_realref(z));

        a = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
        b = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);
        a = fabs(a);
        b = fabs(b);
        absz = sqrt(a * a + b * b);

        if (a <= 1.0 && b <= 1.0)
        {
            use_asymp = 0;
        }
        else if (a > prec || b > prec)
        {
            use_asymp = 1;
        }
        else
        {
            asymp_accuracy = absz * rlog2 * 0.999 - 5;
            use_asymp = asymp_accuracy > prec;
        }

        acb_init(m);
        mag_init(err);

        /* assumes that we already handled the branch cut */
        acb_hypgeom_si_prop_error(err, z);
        acb_get_mid(m, z);

        if (use_asymp)
        {
            acb_hypgeom_si_asymp(res, m, prec);
        }
        else
        {
            /* terms grow to ~ exp(|z|), sum is ~ exp(|im(z)|) */
            cancellation = (absz - b) * rlog2;
            wp = prec + FLINT_MAX(0, cancellation);
            wp = wp * 1.001 + 5;
            acb_hypgeom_si_1f2(res, m, wp);
            acb_set_round(res, res, prec);
        }

        if (pure_imag)
            arb_add_error_mag(acb_imagref(res), err);
        else
            acb_add_error_mag(res, err);

        acb_clear(m);
        mag_clear(err);
    }
}
