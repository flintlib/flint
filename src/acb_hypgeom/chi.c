/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_chi_asymp(acb_t res, const acb_t z, slong prec)
{
    acb_t t, u, v, one;

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(one);

    acb_one(one);

    /* u = U(1,1,z) */
    acb_hypgeom_u_asymp(u, one, one, z, -1, prec);
    /* v = e^(-z) */
    acb_neg(v, z);
    acb_exp(v, v, prec);
    acb_mul(t, u, v, prec);

    if (arb_is_zero(acb_realref(z)))
    {
        arb_div(acb_realref(t), acb_imagref(t), acb_imagref(z), prec);
        arb_zero(acb_imagref(t));
        acb_neg(t, t);
    }
    else
    {
        /* u = U(1,1,-z) */
        acb_neg(u, z);
        acb_hypgeom_u_asymp(u, one, one, u, -1, prec);
        acb_inv(v, v, prec);
        acb_submul(t, u, v, prec);

        acb_div(t, t, z, prec);
        acb_mul_2exp_si(t, t, -1);
        acb_neg(t, t);
    }

    if (acb_is_real(z))
    {
        if (arb_is_positive(acb_realref(z)))
        {
            arb_zero(acb_imagref(t));
        }
        else if (arb_is_negative(acb_realref(z)))
        {
            arb_const_pi(acb_imagref(t), prec);
        }
        else
        {
            /* add [-pi,pi]/2 i */
            acb_const_pi(u, prec);
            arb_zero(acb_imagref(t));
            arb_add_error(acb_imagref(t), acb_realref(u));
        }
    }
    else
    {
        /* -pi/2 if positive real or in lower half plane
           pi/2 if negative real or in upper half plane */
        if (arb_is_negative(acb_imagref(z)))
        {
            acb_const_pi(u, prec);
            acb_mul_2exp_si(u, u, -1);
            arb_sub(acb_imagref(t), acb_imagref(t), acb_realref(u), prec);
        }
        else if (arb_is_positive(acb_imagref(z)))
        {
            acb_const_pi(u, prec);
            acb_mul_2exp_si(u, u, -1);
            arb_add(acb_imagref(t), acb_imagref(t), acb_realref(u), prec);
        }
        else
        {
            /* add [-pi,pi]/2 i */
            acb_const_pi(u, prec);
            acb_mul_2exp_si(u, u, -1);
            arb_add_error(acb_imagref(t), acb_realref(u));
        }
    }

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(one);
}

/* defined in ci.c */
void
acb_hypgeom_chi_2f3(acb_t res, const acb_t z, slong prec);

/* Bound propagated error |chi(z)-chi(mid(z))| assuming that we
   are off the branch cut (not checked here).
   Uses |chi'(z)| = |cosh(z)/z| <= cosh(|re(z)|)/|z|.
*/
static void
acb_hypgeom_chi_prop_error(mag_t err, const acb_t z)
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

        arb_get_mag(t, acb_realref(z));
        mag_cosh(t, t);
        acb_get_mag_lower(u, z);
        mag_div(t, t, u);
        mag_hypot(u, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul(err, t, u);

        mag_clear(t);
        mag_clear(u);
    }
}

void
acb_hypgeom_chi(acb_t res, const acb_t z, slong prec)
{
    if (!acb_is_finite(z) || acb_contains_zero(z))
    {
        acb_indeterminate(res);
    }
/* todo -- good real code
    else if (acb_is_real(z))
    {
    }
*/
    else if (arb_contains_zero(acb_imagref(z)) &&
                !arb_is_nonnegative(acb_imagref(z)) &&
                !arb_is_positive(acb_realref(z)))
    {
        /* We straddle the branch cut; do something generic */
        if (acb_hypgeom_u_use_asymp(z, prec))
            acb_hypgeom_chi_asymp(res, z, prec);
        else
            acb_hypgeom_chi_2f3(res, z, prec);
    }
    else
    {
        double asymp_accuracy, a, b, absz, cancellation, rlog2 = 1.4426950408889634;
        slong wp;
        int use_asymp;
        acb_t m;
        mag_t err;
        /* since we don't have a special case for real z */
        int pure_real = arb_is_zero(acb_imagref(z));

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
        acb_hypgeom_chi_prop_error(err, z);
        acb_get_mid(m, z);

        if (use_asymp)
        {
            acb_hypgeom_chi_asymp(res, m, prec);
        }
        else
        {
            /* terms grow to ~ exp(|z|), sum is ~ exp(|re(z)|) */
            cancellation = (absz - a) * rlog2;
            wp = prec + FLINT_MAX(0, cancellation);
            wp = wp * 1.001 + 5;
            acb_hypgeom_chi_2f3(res, m, wp);
            acb_set_round(res, res, prec);
        }

        if (pure_real)
            arb_add_error_mag(acb_realref(res), err);
        else
            acb_add_error_mag(res, err);

        acb_clear(m);
        mag_clear(err);
    }
}
