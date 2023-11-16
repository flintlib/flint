/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void mag_agm(mag_t res, const mag_t x, const mag_t y);

static void
agm_helper(acb_t res, const acb_t a, const acb_t b, slong prec)
{
    if (acb_rel_accuracy_bits(b) >= acb_rel_accuracy_bits(a))
    {
        acb_div(res, a, b, prec);
        acb_agm1(res, res, prec);
        acb_mul(res, res, b, prec);
    }
    else
    {
        acb_div(res, b, a, prec);
        acb_agm1(res, res, prec);
        acb_mul(res, res, a, prec);
    }
}

void
acb_agm(acb_t res, const acb_t a, const acb_t b, slong prec)
{
    acb_t t, u, v;

    if (!acb_is_finite(a) || !acb_is_finite(b))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(a) || acb_is_zero(b))
    {
        acb_zero(res);
        return;
    }

    if (arb_is_zero(acb_imagref(a)) && arb_is_zero(acb_imagref(b)))
    {
        if (arb_is_nonnegative(acb_realref(a)) && arb_is_nonnegative(acb_realref(b)))
        {
            arb_agm(acb_realref(res), acb_realref(a), acb_realref(b), prec);
            arb_zero(acb_imagref(res));
            return;
        }
    }

    if (acb_contains_zero(a) || acb_contains_zero(b))
    {
        mag_t ra, rb;

        mag_init(ra);
        mag_init(rb);

        acb_get_mag(ra, a);
        acb_get_mag(rb, b);
        mag_agm(ra, ra, rb);
        acb_zero(res);
        acb_add_error_mag(res, ra);

        mag_clear(ra);
        mag_clear(rb);

        return;
    }

    acb_init(t);

    acb_add(t, a, b, prec);
    acb_mul_2exp_si(t, t, -1);

    /* a ~= -b; bound magnitude */
    if (acb_contains_zero(t))
    {
        mag_t ra, rb;

        mag_init(ra);
        mag_init(rb);

        acb_get_mag(ra, a);
        acb_get_mag(rb, b);
        mag_mul(rb, ra, rb);
        mag_sqrt(rb, rb);

        acb_get_mag(ra, t);
        mag_agm(ra, ra, rb);

        acb_zero(res);
        acb_add_error_mag(res, ra);

        mag_clear(ra);
        mag_clear(rb);

        acb_clear(t);
        return;
    }

    /* Do the initial step with the optimal square root, reducing to agm1 */

    acb_init(u);
    acb_init(v);

    acb_mul(u, a, b, prec);

    /* we can compute either square root here; avoid the branch cut */
    if (arf_sgn(arb_midref(acb_realref(u))) >= 0)
    {
        acb_sqrt(u, u, prec);
    }
    else
    {
        acb_neg(u, u);
        acb_sqrt(u, u, prec);
        acb_mul_onei(u, u);
    }

    acb_div(v, t, u, prec);

    if (arb_is_nonnegative(acb_realref(v)))
    {
        agm_helper(res, t, u, prec);
    }
    else if (arb_is_negative(acb_realref(v)))
    {
        acb_neg(u, u);
        agm_helper(res, t, u, prec);
    }
    else
    {
        agm_helper(v, t, u, prec);
        acb_neg(u, u);
        agm_helper(res, t, u, prec);
        acb_union(res, res, v, prec);
    }

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

