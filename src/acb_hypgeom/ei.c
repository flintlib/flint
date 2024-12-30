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
acb_hypgeom_ei_asymp(acb_t res, const acb_t z, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    acb_one(t);
    acb_neg(u, z);

    acb_hypgeom_u_asymp(u, t, t, u, -1, prec);
    acb_div(u, u, z, prec);
    acb_exp(t, z, prec);
    acb_mul(u, u, t, prec);

    if (arb_is_zero(acb_imagref(z)))
    {
        arb_zero(acb_imagref(u));
    }
    else if (arb_is_positive(acb_imagref(z)))
    {
        acb_const_pi(t, prec);
        arb_add(acb_imagref(u), acb_imagref(u), acb_realref(t), prec);
    }
    else if (arb_is_negative(acb_imagref(z)))
    {
        acb_const_pi(t, prec);
        arb_sub(acb_imagref(u), acb_imagref(u), acb_realref(t), prec);
    }
    else
    {
        /* add [-pi,pi] i */
        acb_const_pi(t, prec);
        arb_add_error(acb_imagref(u), acb_realref(t));
    }

    acb_swap(res, u);

    acb_clear(t);
    acb_clear(u);
}

/*
Ei(z) = z 2F2(1,1,2,2,z) + 0.5[log(z)-log(1/z)] + gamma
*/
void
acb_hypgeom_ei_2f2(acb_t res, const acb_t z, slong prec)
{
    acb_t a, t;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(t);

    acb_one(a);
    acb_set_ui(b, 2);
    acb_set_ui(b + 1, 2);

    acb_hypgeom_pfq_direct(t, a, 1, b, 2, z, -1, prec);
    acb_mul(t, t, z, prec);

    arb_const_euler(acb_realref(a), prec);
    arb_add(acb_realref(t), acb_realref(t), acb_realref(a), prec);

    if (arb_is_zero(acb_imagref(z)))
    {
        if (arb_is_positive(acb_realref(z)))
        {
            acb_log(a, z, prec);
        }
        else
        {
            /* ok if overlapping zero -- will be indeterminate either way */
            acb_neg(a, z);
            acb_log(a, a, prec);
            arb_zero(acb_imagref(a));
        }

        acb_add(t, t, a, prec);
    }
    else if (arb_is_nonzero(acb_imagref(z))
        || arb_is_nonnegative(acb_realref(z))) /* not overlapping (-inf,0] */
    {
        acb_log(a, z, prec);
        acb_add(t, t, a, prec);
    }
    else
    {
        acb_log(a, z, prec);
        arb_zero(acb_imagref(a));

        acb_const_pi(b, prec);
        arb_add_error(acb_imagref(a), acb_realref(b));

        acb_add(t, t, a, prec);
    }

    acb_swap(res, t);

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(t);
}

/* Bound propagated error |ei(z)-ei(mid(z))| assuming that we
   are off the branch cut (not checked here).
   Uses |ei'(z)| = |exp(z)/z| <= exp(re(z))/|z|.
*/
static void
acb_hypgeom_ei_prop_error(mag_t err, const acb_t z)
{
    if (acb_is_exact(z))
    {
        mag_zero(err);
    }
    else
    {
        mag_t t, u;
        arf_t x;

        mag_init(t);
        mag_init(u);
        arf_init(x);

        arb_get_ubound_arf(x, acb_realref(z), MAG_BITS);

        if (arf_sgn(x) >= 0)
        {
            arf_get_mag(t, x);
            mag_exp(t, t);
        }
        else
        {
            arf_get_mag_lower(t, x);
            mag_expinv(t, t);
        }

        acb_get_mag_lower(u, z);
        mag_div(t, t, u);
        mag_hypot(u, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul(err, t, u);

        mag_clear(t);
        mag_clear(u);
        arf_clear(x);
    }
}

void
acb_hypgeom_ei(acb_t res, const acb_t z, slong prec)
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
                !arb_is_zero(acb_imagref(z)) &&
                !arb_is_positive(acb_realref(z)))
    {
        /* We straddle the branch cut; do something generic */
        if (acb_hypgeom_u_use_asymp(z, prec))
            acb_hypgeom_ei_asymp(res, z, prec);
        else
            acb_hypgeom_ei_2f2(res, z, prec);
    }
    else
    {
        double asymp_accuracy, a, asigned, b, absz, cancellation, rlog2 = 1.4426950408889634;
        slong wp;
        int use_asymp;
        acb_t m;
        mag_t err;
        /* since we don't have a special case for real z */
        int pure_real = arb_is_zero(acb_imagref(z));

        a = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
        b = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);
        asigned = a;
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
        acb_hypgeom_ei_prop_error(err, z);
        acb_get_mid(m, z);

        if (use_asymp)
        {
            acb_hypgeom_ei_asymp(res, m, prec);
        }
        else
        {
            /* terms grow to ~ exp(|z|), sum is ~ exp(re(z)) */
            cancellation = (absz - asigned) * rlog2;
            wp = prec + FLINT_MAX(0, cancellation);
            wp = wp * 1.001 + 5;
            acb_hypgeom_ei_2f2(res, m, wp);
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
