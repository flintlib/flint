/*
    Copyright (C) 2013, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "fixed.h"

static void
arb_const_log2_hypgeom_eval(arb_t s, slong prec)
{
    /* Zuniga's series, 11.9 bits per term (fixed/const_log2.c) */
    _fixed_constant_arb(s, fball_const_log2, prec);
}

_ARB_DEF_CACHED_CONSTANT(static, arb_const_log2_hypgeom, arb_const_log2_hypgeom_eval)

void
arb_const_log2(arb_t res, slong prec)
{
    if (prec < ARB_LOG_TAB2_LIMBS * FLINT_BITS - 16)
    {
        slong exp;

        /* just reading the table is known to give the correct rounding */
        _arf_set_round_mpn(arb_midref(res), &exp, arb_log_log2_tab,
            ARB_LOG_TAB2_LIMBS, 0, prec, ARF_RND_NEAR);
        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp);

        /* 1/2 ulp error */
        _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), exp - prec);
        MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
    }
    else
    {
        arb_const_log2_hypgeom(res, prec);
    }
}
