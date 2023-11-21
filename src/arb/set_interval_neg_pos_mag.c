/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

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
            arf_neg(arb_midref(res), arb_midref(res));
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

