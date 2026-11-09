/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

static void
arb_const_pi_chudnovsky_eval(arb_t s, slong prec)
{
    /* the Chudnovsky series (mp_real/const_pi.c) */
    _mp_real_const_arb(s, mp_real_const_pi4, prec);
    arb_mul_2exp_si(s, s, 2);
}

_ARB_DEF_CACHED_CONSTANT(static, arb_const_pi_chudnovsky, arb_const_pi_chudnovsky_eval)

void
arb_const_pi(arb_t res, slong prec)
{
    if (prec <= FLINT_BITS * (MP_REAL_CONST_STATIC_N - 1) - 16)
        _arb_const_set_floor(res, MP_REAL_CONST_ID_PI4, 2, prec);
    else
        _arb_const_round_cached(res, arb_const_pi_chudnovsky, prec);
}
