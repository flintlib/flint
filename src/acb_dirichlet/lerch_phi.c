/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"
#include "acb_dirichlet.h"

void
acb_dirichlet_lerch_phi(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
{
    if (!acb_is_finite(z) || !acb_is_finite(s) || !acb_is_finite(a))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_contains_int(a) && !arb_is_positive(acb_realref(a)))
    {
        if (!(acb_is_int(s) && arb_is_nonpositive(acb_realref(s))))
        {
            acb_indeterminate(res);
            return;
        }
    }

    if (acb_is_zero(z))
    {
        acb_t t;
        acb_init(t);
        acb_neg(t, s);
        acb_pow(res, a, t, prec);
        acb_clear(t);
        return;
    }

    if (acb_is_one(z))
    {
        arb_t one;
        arb_init(one);
        if (arb_gt(acb_realref(s), one))
            acb_dirichlet_hurwitz(res, s, a, prec);
        else
            acb_indeterminate(res);
        arb_clear(one);
        return;
    }

    if (acb_equal_si(z, -1))
    {
        if (acb_is_one(a))
        {
            acb_dirichlet_eta(res, s, prec);
        }
        else if (acb_is_one(s))
        {
            /* (psi((a+1)/2) - psi(a/2))/2 */
            acb_t t, u;
            acb_init(t);
            acb_init(u);
            acb_mul_2exp_si(t, a, -1);
            acb_digamma(t, t, prec);
            acb_add_ui(u, a, 1, prec);
            acb_mul_2exp_si(u, u, -1);
            acb_digamma(u, u, prec);
            acb_sub(res, u, t, prec);
            acb_mul_2exp_si(res, res, -1);
            acb_clear(t);
            acb_clear(u);
        }
        else
        {
            /* 2^(-s) (zeta(s,a/2) - zeta(s,(a+1)/2)) */
            acb_t t, u;
            acb_init(t);
            acb_init(u);
            acb_mul_2exp_si(t, a, -1);
            acb_hurwitz_zeta(t, s, t, prec);
            acb_add_ui(u, a, 1, prec);
            acb_mul_2exp_si(u, u, -1);
            acb_hurwitz_zeta(u, s, u, prec);
            acb_sub(t, t, u, prec);
            acb_neg(u, s);
            acb_set_ui(res, 2);
            acb_pow(res, res, u, prec);
            acb_mul(res, res, t, prec);
            acb_clear(t);
            acb_clear(u);
        }
        return;
    }

    if (acb_is_zero(s))
    {
        acb_sub_ui(res, z, 1, prec + 5);
        acb_neg(res, res);
        acb_inv(res, res, prec);
        return;
    }

    if (acb_is_one(s))
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);
        acb_one(t);
        acb_add_ui(u, a, 1, prec + 5);
        acb_hypgeom_2f1(t, t, a, u, z, ACB_HYPGEOM_2F1_BC, prec + 5);
        acb_div(res, t, a, prec);
        if (!acb_is_finite(res))
            acb_indeterminate(res);
        acb_clear(t);
        acb_clear(u);
        return;
    }

    if (acb_is_one(a) && !acb_contains_zero(z))
    {
        acb_t t;
        acb_init(t);
        acb_polylog(t, s, z, prec);
        acb_div(res, t, z, prec);
        acb_clear(t);
        return;
    }

    {
        mag_t zm, lim;
        mag_init(zm);
        mag_init(lim);

        acb_get_mag(zm, z);
        mag_set_d(lim, 0.875);

        if (mag_cmp(zm, lim) <= 0)
        {
            acb_dirichlet_lerch_phi_direct(res, z, s, a, prec);
        }
        else
        {
            acb_dirichlet_lerch_phi_integral(res, z, s, a, prec);
        }

        mag_clear(zm);
        mag_clear(lim);
    }
}
