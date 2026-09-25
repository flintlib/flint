/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "fixed.h"

static void
arb_const_pi_chudnovsky_eval(arb_t s, slong prec)
{
    /* the Chudnovsky series (fixed/const_pi.c) */
    _fixed_constant_arb(s, fball_const_pi_chudnovsky, prec);
}

_ARB_DEF_CACHED_CONSTANT(static, arb_const_pi_chudnovsky, arb_const_pi_chudnovsky_eval)

void
arb_const_pi(arb_t res, slong prec)
{
    if (prec < ARB_PI4_TAB_LIMBS * FLINT_BITS - 16)
    {
        slong exp;

        /* just reading the table is known to give the correct rounding */
        _arf_set_round_mpn(arb_midref(res), &exp, arb_pi4_tab,
            ARB_PI4_TAB_LIMBS, 0, prec, ARF_RND_NEAR);
        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), 2 + exp);

        /* 1/2 ulp error */
        _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), 2 + exp - prec);
        MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
    }
    else
    {
        arb_const_pi_chudnovsky(res, prec);
    }
}
