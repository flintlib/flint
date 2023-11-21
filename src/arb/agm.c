/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

/* Use series expansion of the elliptic integral
   pi/(4K(z^2)) = 1/2 - z^2/8 - 5z^4/128 - 11z^6/512 - 469z^8/32768 + O(z^10)
   to avoid computing the last couple of AGM iterations.
   The higher coefficients are bounded by 1/64, so the tail is easily
   bounded by a geometric series (the error is passed as input here).
   Note: the input data is a+b, a-b; the variables z, z2 are scratch space. */
static void
arb_agm_close_taylor(arb_t res, arb_t z, arb_t z2,
    const arb_t aplusb, const arb_t aminusb,
    const mag_t err, slong prec)
{
    arb_div(z, aminusb, aplusb, prec);
    arb_sqr(z, z, prec);
    arb_sqr(z2, z, prec);

    arb_mul_si(res, z2, -469, prec);
    arb_addmul_si(res, z, -704, prec);
    arb_mul(res, res, z2, prec);
    arb_addmul_si(res, z2, -1280, prec);
    arb_mul_2exp_si(z, z, 12);
    arb_sub(res, res, z, prec);
    arb_add_ui(res, res, 16384, prec);
    arb_mul_2exp_si(res, res, -15);

    arb_add_error_mag(res, err);
    arb_mul(res, res, aplusb, prec);
}

void
mag_agm(mag_t res, const mag_t x, const mag_t y)
{
    if (!mag_is_finite(x) || !mag_is_finite(y))
    {
        mag_inf(res);
    }
    else if (mag_is_zero(x) || mag_is_zero(y))
    {
        mag_zero(res);
    }
    else
    {
        mag_t t, u, a, b, one_eps;

        mag_init(a);
        mag_init(b);
        mag_init(t);
        mag_init(u);
        mag_init(one_eps);

        /* invariant: a is an upper bound for agm(x,y) */
        /* b would be a lower bound if not for rounding errors */
        mag_max(a, x, y);
        mag_min(b, x, y);

        mag_one(one_eps);
        mag_add_ui_2exp_si(one_eps, one_eps, 1, -26);

        while (1)
        {
            mag_mul(t, b, one_eps);
            if (mag_cmp(t, a) > 0)
            {
                mag_set(res, a);
                break;
            }

            mag_add(t, a, b);
            mag_mul_2exp_si(t, t, -1);
            mag_mul(u, a, b);
            mag_sqrt(u, u);
            mag_swap(t, a);
            mag_swap(u, b);
        }

        mag_clear(a);
        mag_clear(b);
        mag_clear(t);
        mag_clear(u);
        mag_clear(one_eps);
    }
}

void
mag_agm_lower(mag_t res, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x) || mag_is_zero(y))
    {
        mag_zero(res);
    }
    else if (!mag_is_finite(x) || !mag_is_finite(y))
    {
        mag_inf(res);
    }
    else
    {
        mag_t t, u, a, b, one_eps;

        mag_init(a);
        mag_init(b);
        mag_init(t);
        mag_init(u);
        mag_init(one_eps);

        /* invariant: b is a lower bound for agm(x,y) */
        /* a would be an upper bound if not for rounding errors */
        mag_max(a, x, y);
        mag_min(b, x, y);

        mag_one(one_eps);
        mag_add_ui_2exp_si(one_eps, one_eps, 1, -26);

        while (1)
        {
            mag_mul(t, b, one_eps);
            if (mag_cmp(t, a) > 0)
            {
                mag_set(res, b);
                break;
            }

            mag_add_lower(t, a, b);
            mag_mul_2exp_si(t, t, -1);
            mag_mul_lower(u, a, b);
            mag_sqrt_lower(u, u);
            mag_swap(t, a);
            mag_swap(u, b);
        }

        mag_clear(a);
        mag_clear(b);
        mag_clear(t);
        mag_clear(u);
        mag_clear(one_eps);
    }
}


void
arb_agm(arb_t res, const arb_t x, const arb_t y, slong prec)
{
    arb_t a, b, t, u;
    mag_t err, err2;
    slong acc1, acc2, wp;

    if (!arb_is_nonnegative(x) || !arb_is_nonnegative(y) ||
        !arb_is_finite(x) || !arb_is_finite(y))
    {
        arb_indeterminate(res);
        return;
    }

    if (arb_is_zero(x) || arb_is_zero(y))
    {
        arb_zero(res);
        return;
    }

    arb_init(a);
    arb_init(b);
    arb_init(t);
    arb_init(u);
    mag_init(err);
    mag_init(err2);

    arb_set(a, x);
    arb_set(b, y);

    wp = prec;

    while (1)
    {
        acc1 = arb_rel_accuracy_bits(a);
        acc2 = arb_rel_accuracy_bits(b);
        acc1 = FLINT_MIN(acc1, acc2);

        /* Compute lower and upper bounds if we don't need high precision. */
        if (acc1 < 20 || wp < 20)
        {
            arb_get_mag_lower(arb_radref(t), a);
            arb_get_mag_lower(arb_radref(u), b);
            mag_agm_lower(err, arb_radref(t), arb_radref(u));

            arb_get_mag(arb_radref(t), a);
            arb_get_mag(arb_radref(u), b);
            mag_agm(err2, arb_radref(t), arb_radref(u));

            arf_set_mag(arb_midref(t), err);
            arf_set_mag(arb_midref(u), err2);

            arb_set_interval_arf(res, arb_midref(t), arb_midref(u), prec);
            break;
        }

        if (acc1 < wp - 2 * MAG_BITS)
            wp = acc1 + 2 * MAG_BITS;

        arb_sub(u, a, b, wp);

        /* Fallback exit. */
        if (arb_contains_zero(u))
        {
            arb_union(res, a, b, wp);
            break;
        }

        arb_add(t, a, b, wp);

        arb_get_mag(err, u);
        arb_get_mag_lower(err2, t);
        mag_div(err, err, err2);
        mag_geom_series(err, err, 10);
        mag_mul_2exp_si(err, err, -6);

        /* Use Taylor series when we have 1/10 the accurate bits. */
        if (mag_cmp_2exp_si(err, -wp) < 0)
        {
            /* pass a, b as scratch space */
            arb_agm_close_taylor(res, a, b, t, u, err, wp);
            break;
        }

        arb_mul_2exp_si(t, t, -1);
        arb_mul(u, a, b, wp);
        arb_sqrt(u, u, wp);
        arb_swap(t, a);
        arb_swap(u, b);
    }

    arb_clear(a);
    arb_clear(b);
    arb_clear(t);
    arb_clear(u);
    mag_clear(err);
    mag_clear(err2);
}

