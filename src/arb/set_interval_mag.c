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

