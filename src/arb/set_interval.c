/*
    Copyright (C) 2013, 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arb.h"

void
arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, slong prec)
{
    arf_t t;
    int inexact;

    /* [-inf, -inf] or [+inf, +inf] */
    if (arf_is_inf(a) && arf_equal(a, b))
    {
        arf_set(arb_midref(x), a);
        mag_zero(arb_radref(x));
        return;
    }

    /* any nan -> [nan +/- inf] */
    if (arf_is_nan(a) || arf_is_nan(b))
    {
        arb_indeterminate(x);
        return;
    }

    /* [-inf, x] or [x, +inf] = [+/- inf] */
    if (arf_is_neg_inf(a) || arf_is_pos_inf(b))
    {
        arf_zero(arb_midref(x));
        mag_inf(arb_radref(x));
        return;
    }

    /* [(a + b) +/- (b - a)] / 2 */
    arf_init(t);
    arf_sub(t, b, a, MAG_BITS, ARF_RND_UP);

    if (arf_sgn(t) < 0)
    {
        flint_throw(FLINT_ERROR, "exception: arb_set_interval_arf: endpoints not ordered\n");
    }

    arf_get_mag(arb_radref(x), t);

    inexact = arf_add(arb_midref(x), a, b, prec, ARB_RND);
    if (inexact)
        arf_mag_add_ulp(arb_radref(x), arb_radref(x), arb_midref(x), prec);

    arb_mul_2exp_si(x, x, -1);

    arf_clear(t);
}

void
arb_set_interval_mag(arb_t res, const mag_t a, const mag_t b, slong prec)
{
    if (MAG_IS_LAGOM(a) && MAG_IS_LAGOM(b))
    {
        slong aexp, bexp;
        mp_limb_t aman, bman, mman, rman, tmp;

        aman = MAG_MAN(a);
        bman = MAG_MAN(b);

        aexp = MAG_EXP(a);
        bexp = MAG_EXP(b);

        if (aman == 0 && bman == 0)
        {
            arb_zero(res);
            return;
        }

        if (bman == 0 || (aman != 0 &&
                            (aexp > bexp || (aexp == bexp && aman > bman))))
        {
            flint_throw(FLINT_ERROR, "exception: arb_set_interval_mag: endpoints not ordered\n");
        }

        /* now a = 0 or bexp >= aexp */
        if (aman == 0 || bexp - aexp > MAG_BITS)
        {
            mman = bman;                     /* midpoint a+b */
            rman = bman + (aman != 0);       /* radius b-a */
        }
        else
        {
            tmp = (aman >> (bexp - aexp));
            mman = bman + tmp;                         /* midpoint a+b */
            rman = bman - tmp;                         /* radius b-a */
            rman += ((tmp << (bexp - aexp)) != aman);  /* rounding error */
        }

        arf_set_ui(arb_midref(res), mman);
        /* m can't be zero */
        ARF_EXP(arb_midref(res)) += bexp - MAG_BITS - 1;

        mag_set_ui(arb_radref(res), rman);
        if (rman != 0)  /* r can be zero */
            MAG_EXP(arb_radref(res)) += bexp - MAG_BITS - 1;

        arb_set_round(res, res, prec);
    }
    else
    {
        int inexact;
        arf_t aa, bb;

        if (mag_cmp(a, b) > 0)
        {
            flint_throw(FLINT_ERROR, "exception: arb_set_interval_mag: endpoints not ordered\n");
        }

        if (mag_is_inf(a))
        {
            arb_pos_inf(res);
            return;
        }

        if (mag_is_inf(b))
        {
            arb_zero_pm_inf(res);
            return;
        }

        arf_init_set_mag_shallow(aa, a);
        arf_init_set_mag_shallow(bb, b);

        inexact = arf_add(arb_midref(res), aa, bb, prec, ARB_RND);

        mag_sub(arb_radref(res), b, a);

        if (inexact)
            arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);

        arb_mul_2exp_si(res, res, -1);
    }
}

void
arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, slong prec)
{
    arf_t aa, bb;

    arf_init(aa);
    arf_init(bb);

    arf_set_mpfr(aa, a);
    arf_set_mpfr(bb, b);

    arb_set_interval_arf(x, aa, bb, prec);

    arf_clear(aa);
    arf_clear(bb);
}

void
arb_set_interval_neg_pos_mag(arb_t res, const mag_t a, const mag_t b, slong prec)
{
    if (MAG_IS_LAGOM(a) && MAG_IS_LAGOM(b))
    {
        slong aexp, bexp, mexp, shift;
        mp_limb_t aman, bman, mman, rman, tmp;
        int negative;

        aman = MAG_MAN(a);
        bman = MAG_MAN(b);

        aexp = MAG_EXP(a);
        bexp = MAG_EXP(b);

        if (aman == 0)
        {
            if (bman == 0)
            {
                arb_zero(res);
                return;
            }

            negative = 0;
            mexp = bexp;
            mman = bman;
            rman = bman;
        }
        else if (bman == 0)
        {
            negative = 1;
            mexp = aexp;
            mman = aman;
            rman = aman;
        }
        else if (aexp == bexp)
        {
            mexp = aexp;
            negative = aman >= bman;
            if (negative)
                mman = aman - bman;
            else
                mman = bman - aman;
            rman = aman + bman;
        }
        else if (aexp > bexp)
        {
            negative = 1;
            mexp = aexp;
            shift = aexp - bexp;
            if (shift > MAG_BITS)
            {
                mman = aman;
                rman = aman + 2;
            }
            else
            {
                tmp = bman >> shift;
                mman = aman - tmp;
                rman = aman + tmp;
                rman += 2 * ((tmp << shift) != bman);
            }
        }
        else
        {
            negative = 0;
            mexp = bexp;
            shift = bexp - aexp;
            if (shift > MAG_BITS)
            {
                mman = bman;
                rman = bman + 2;
            }
            else
            {
                tmp = aman >> shift;
                mman = bman - tmp;
                rman = bman + tmp;
                rman += 2 * ((tmp << shift) != aman);
            }
        }

        arf_set_ui(arb_midref(res), mman);
        if (negative)
            arf_inplace_neg(arb_midref(res));
        if (mman != 0)
            ARF_EXP(arb_midref(res)) += mexp - MAG_BITS - 1;

        mag_set_ui(arb_radref(res), rman);
        /* r can't be zero */
        MAG_EXP(arb_radref(res)) += mexp - MAG_BITS - 1;

        arb_set_round(res, res, prec);
    }
    else
    {
        arf_t aa, bb;
        int inexact;

        if (mag_is_inf(a) || mag_is_inf(b))
        {
            arb_zero_pm_inf(res);
            return;
        }

        arf_init_set_mag_shallow(aa, a);
        arf_init_set_mag_shallow(bb, b);

        inexact = arf_sub(arb_midref(res), bb, aa, prec, ARB_RND);

        mag_add(arb_radref(res), b, a);

        if (inexact)
            arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec);

        arb_mul_2exp_si(res, res, -1);
    }
}
