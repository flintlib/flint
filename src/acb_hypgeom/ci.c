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
acb_hypgeom_ci_asymp(acb_t res, const acb_t z, slong prec)
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
        arb_div(acb_realref(t), acb_imagref(t), acb_realref(z), prec);
        arb_zero(acb_imagref(t));
        acb_neg(t, t);
    }
    else
    {
        /* u = U(1,1,-iz) */
        acb_neg(w, w);
        acb_hypgeom_u_asymp(u, one, one, w, -1, prec);
        acb_inv(v, v, prec);
        acb_submul(t, u, v, prec);

        acb_div(t, t, w, prec);
        acb_mul_2exp_si(t, t, -1);
    }

    if (arb_is_zero(acb_realref(z)))
    {
        if (arb_is_positive(acb_imagref(z)))
        {
            arb_const_pi(acb_imagref(t), prec);
            arb_mul_2exp_si(acb_imagref(t), acb_imagref(t), -1);
        }
        else if (arb_is_negative(acb_imagref(z)))
        {
            arb_const_pi(acb_imagref(t), prec);
            arb_mul_2exp_si(acb_imagref(t), acb_imagref(t), -1);
            arb_neg(acb_imagref(t), acb_imagref(t));
        }
        else
        {
            acb_const_pi(u, prec);
            acb_mul_2exp_si(u, u, -1);
            arb_zero(acb_imagref(t));
            arb_add_error(acb_imagref(t), acb_realref(u));
        }
    }
    else
    {
        /* 0 if positive or positive imaginary
           pi if upper left quadrant (including negative real axis)
           -pi if lower left quadrant (including negative imaginary axis) */
        if (arb_is_positive(acb_realref(z)))
        {
            /* do nothing */
        }
        else if (arb_is_negative(acb_realref(z)) && arb_is_nonnegative(acb_imagref(z)))
        {
            acb_const_pi(u, prec);
            arb_add(acb_imagref(t), acb_imagref(t), acb_realref(u), prec);
        }
        else if (arb_is_nonpositive(acb_realref(z)) && arb_is_negative(acb_imagref(z)))
        {
            acb_const_pi(u, prec);
            arb_sub(acb_imagref(t), acb_imagref(t), acb_realref(u), prec);
        }
        else
        {
            /* add [-pi,pi] */
            acb_const_pi(u, prec);
            arb_add_error(acb_imagref(t), acb_realref(u));
        }
    }

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(w);
    acb_clear(v);
    acb_clear(one);
}

static void
_acb_hypgeom_ci_2f3(acb_t res, const acb_t z, int hyperbolic, slong prec)
{
    acb_t a, t, u;
    acb_struct b[3];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(b + 2);
    acb_init(t);
    acb_init(u);

    acb_one(a);
    acb_set_ui(b, 2);
    acb_set(b + 1, b);
    acb_set_ui(b + 2, 3);
    acb_mul_2exp_si(b + 2, b + 2, -1);

    acb_mul(t, z, z, prec);
    acb_mul_2exp_si(t, t, -2);
    if (!hyperbolic)
        acb_neg(t, t);
    acb_hypgeom_pfq_direct(u, a, 1, b, 3, t, -1, prec);
    acb_mul(u, u, t, prec);

    acb_log(t, z, prec);
    acb_add(u, u, t, prec);

    arb_const_euler(acb_realref(t), prec);
    arb_add(acb_realref(u), acb_realref(u), acb_realref(t), prec);

    acb_swap(res, u);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(b + 2);
    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_ci_2f3(acb_t res, const acb_t z, slong prec)
{
    _acb_hypgeom_ci_2f3(res, z, 0, prec);
}

void
acb_hypgeom_chi_2f3(acb_t res, const acb_t z, slong prec)
{
    _acb_hypgeom_ci_2f3(res, z, 1, prec);
}

/* Bound propagated error |ci(z)-ci(mid(z))| assuming that we
   are off the branch cut (not checked here).
   Uses |ci'(z)| = |cos(z)/z| <= cosh(|im(z)|)/|z|.
*/
static void
acb_hypgeom_ci_prop_error(mag_t err, const acb_t z)
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
        mag_div(t, t, u);
        mag_hypot(u, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul(err, t, u);

        mag_clear(t);
        mag_clear(u);
    }
}

void
acb_hypgeom_ci(acb_t res, const acb_t z, slong prec)
{
    if (!acb_is_finite(z) || acb_contains_zero(z))
    {
        acb_indeterminate(res);
    }
    else if (acb_is_real(z))
    {
        if (arb_is_positive(acb_realref(z)))
        {
            arb_hypgeom_ci(acb_realref(res), acb_realref(z), prec);
            arb_zero(acb_imagref(res));
        }
        else /* arb_is_negative(acb_realref(z)), since we already tested for 0 */
        {
            arb_neg(acb_realref(res), acb_realref(z));
            arb_hypgeom_ci(acb_realref(res), acb_realref(res), prec);
            arb_const_pi(acb_imagref(res), prec);
        }
    }
    else if (arb_contains_zero(acb_imagref(z)) &&
                !arb_is_nonnegative(acb_imagref(z)) &&
                !arb_is_positive(acb_realref(z)))
    {
        /* We straddle the branch cut; do something generic */
        if (acb_hypgeom_u_use_asymp(z, prec))
            acb_hypgeom_ci_asymp(res, z, prec);
        else
            acb_hypgeom_ci_2f3(res, z, prec);
    }
    else
    {
        double asymp_accuracy, a, b, absz, cancellation, rlog2 = 1.4426950408889634;
        slong wp;
        int use_asymp;
        acb_t m;
        mag_t err;

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
        acb_hypgeom_ci_prop_error(err, z);
        acb_get_mid(m, z);

        if (use_asymp)
        {
            acb_hypgeom_ci_asymp(res, m, prec);
        }
        else
        {
            /* terms grow to ~ exp(|z|), sum is ~ exp(|im(z)|) */
            cancellation = (absz - b) * rlog2;
            wp = prec + FLINT_MAX(0, cancellation);
            wp = wp * 1.001 + 5;
            acb_hypgeom_ci_2f3(res, m, wp);
            acb_set_round(res, res, prec);
        }

        acb_add_error_mag(res, err);

        acb_clear(m);
        mag_clear(err);
    }
}
