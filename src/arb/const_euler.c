/*
    Copyright (C) 2012, 2013, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_hypgeom.h"
#include "mp_real.h"

static void
arb_const_euler_eval(arb_t res, slong prec)
{
    /* Brent-McMillan in dual numbers (mp_real/const_euler.c) */
    _mp_real_const_arb(res, mp_real_const_euler, prec);
}

_ARB_DEF_CACHED_CONSTANT(static, arb_const_euler_brent_mcmillan, arb_const_euler_eval)

FLINT_DLL extern const ulong arb_hypgeom_gamma_tab_limbs[];

void
arb_const_euler(arb_t res, slong prec)
{
    if (prec < ARB_HYPGEOM_GAMMA_TAB_PREC - 16)
    {
        slong exp;
        slong n;

        n = ARB_HYPGEOM_GAMMA_TAB_PREC / FLINT_BITS;

        /* just reading the table is known to give the correct rounding */
        _arf_set_round_mpn(arb_midref(res), &exp, arb_hypgeom_gamma_tab_limbs + n, n, 0, prec, ARF_RND_NEAR);
        _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp);

        /* 1/2 ulp error */
        _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), exp - prec);
        MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
    }
    else
    {
        arb_const_euler_brent_mcmillan(res, prec);
    }
}
