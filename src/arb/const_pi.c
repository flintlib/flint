/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "arb.h"
#include "hypgeom.h"

void
arb_const_pi_chudnovsky_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t, u;

    arb_init(t);
    arb_init(u);
    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "2  13591409 545140134");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "4  5 -46 108 -72");
    fmpz_poly_set_str(series->Q, "4  0 0 0 10939058860032000");

    prec += FLINT_CLOG2(prec) + 5;
    arb_hypgeom_infsum(s, t, series, prec, prec);

    arb_rsqrt_ui(u, 640320, prec);
    arb_mul(s, s, u, prec);

    arb_mul_ui(t, t, 640320 / 12, prec);
    arb_div(s, t, s, prec);

    hypgeom_clear(series);
    arb_clear(t);
    arb_clear(u);
}

ARB_DEF_CACHED_CONSTANT(arb_const_pi_chudnovsky, arb_const_pi_chudnovsky_eval)

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
