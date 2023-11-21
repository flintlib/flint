/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

static void
arb_log_abs(arb_t res, const arb_t a, slong prec)
{
    if (arb_is_positive(a))
    {
        arb_log(res, a, prec);
    }
    else if (arb_is_negative(a))
    {
        arb_neg(res, a);
        arb_log(res, res, prec);
    }
    else
    {
        arb_indeterminate(res);
    }
}

static int
arf_close_to_one(const arf_t z)
{
    mp_limb_t top;

    if (ARF_EXP(z) == 0)
    {
        ARF_GET_TOP_LIMB(top, z);
        return (top >> (FLINT_BITS - 4)) == 15;
    }
    else if (ARF_EXP(z) == 1)
    {
        ARF_GET_TOP_LIMB(top, z);
        return (top >> (FLINT_BITS - 4)) == 8;
    }

    return 0;
}

void
arb_log_hypot(arb_t res, const arb_t a, const arb_t b, slong prec)
{
    slong acc;
    arb_t x;

    if (arb_is_zero(b))
    {
        arb_log_abs(res, a, prec);
        return;
    }

    if (arb_is_zero(a))
    {
        arb_log_abs(res, b, prec);
        return;
    }

    if (!arb_is_finite(a) || !arb_is_finite(b))
    {
        if (arf_is_nan(arb_midref(a)) || arf_is_nan(arb_midref(b)))
        {
            arb_indeterminate(res);
        }
        else if ((!arb_is_finite(a) && !arb_contains_zero(a)) ||
                 (!arb_is_finite(b) && !arb_contains_zero(b)))
        {
            arb_pos_inf(res);
        }
        else
        {
            arb_indeterminate(res);
        }
        return;
    }

    /* a close to 1 -- for accurate arb_log1p */
    if (mag_cmp_2exp_si(arb_radref(a), -3) < 0 &&
        mag_cmp_2exp_si(arb_radref(b), -3) < 0 &&
        arf_cmpabs_2exp_si(arb_midref(b), -3) < 0 &&
        arf_close_to_one(arb_midref(a)))
    {
        arb_t y;
        arb_init(x);
        arb_init(y);
        if (arf_sgn(arb_midref(a)) > 0)
        {
            arb_sub_ui(y, a, 1, prec + 8);
            arb_mul(x, y, y, prec + 8);
            arb_addmul(x, b, b, prec + 8);
            arb_mul_2exp_si(y, y, 1);
            arb_add(x, x, y, prec + 8);
        }
        else
        {
            arb_add_ui(y, a, 1, prec + 8);
            arb_mul(x, y, y, prec + 8);
            arb_addmul(x, b, b, prec + 8);
            arb_mul_2exp_si(y, y, 1);
            arb_sub(x, x, y, prec + 8);
        }
        arb_log1p(res, x, prec);
        arb_mul_2exp_si(res, res, -1);
        arb_clear(x);
        arb_clear(y);
        return;
    }

    arb_init(x);

    /* todo: write an arb_sosq function */
    /* todo: for very wide input, we could predict that a^2+b^2 will have low
       accuracy without computing it and go more quickly to the interval case
       -- however, a first failed attempt to write such code proved that
          it's actually somewhat complicated to do... */
    arb_mul(x, a, a, FLINT_MAX(MAG_BITS, prec) + 8);
    arb_addmul(x, b, b, FLINT_MAX(MAG_BITS, prec) + 8);

    acc = arb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    acc = FLINT_MIN(acc, prec);

    if (acc > 10)
    {
        arb_log(res, x, prec);
        arb_mul_2exp_si(res, res, -1);
    }
    else
    {
        mag_t t, u, v;

        mag_init(t);
        mag_init(u);
        mag_init(v);

        arb_get_mag_lower(t, a);
        arb_get_mag_lower(u, b);

        if (!arb_contains_zero(x))
            acc += fmpz_bits(ARF_EXPREF(arb_midref(x)));

        if (mag_is_zero(t) && mag_is_zero(u))
        {
            arb_indeterminate(res);
        }
        else if (acc < 20)
        {
            /* t = lower bound for a^2 + b^2 */
            mag_mul_lower(t, t, t);
            mag_mul_lower(u, u, u);
            mag_add_lower(t, t, u);

            /* u = upper bound for a^2 + b^2 */
            arb_get_mag(u, x);

            if (mag_cmp_2exp_si(t, 0) >= 0)
            {
                mag_log_lower(t, t);
                mag_log(u, u);
                arb_set_interval_mag(res, t, u, prec);
            }
            else if (mag_cmp_2exp_si(u, 0) <= 0)
            {
                mag_neg_log_lower(u, u);
                mag_neg_log(t, t);
                arb_set_interval_mag(res, u, t, prec);
                arb_neg(res, res);
            }
            else
            {
                mag_neg_log(t, t);
                mag_log(u, u);
                arb_set_interval_neg_pos_mag(res, t, u, prec);
            }

            arb_mul_2exp_si(res, res, -1);
        }
        else
        {
            arb_log(res, x, prec);
            arb_mul_2exp_si(res, res, -1);
        }

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
    }

    arb_clear(x);
}

