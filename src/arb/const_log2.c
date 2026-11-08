/*
    Copyright (C) 2013, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

static void
arb_const_log2_hypgeom_eval(arb_t s, slong prec)
{
    /* Zuniga's series, 11.9 bits per term (fixed/const_log2.c) */
    _mp_real_const_arb(s, mp_real_const_log2, prec);
}

_ARB_DEF_CACHED_CONSTANT(static, arb_const_log2_hypgeom, arb_const_log2_hypgeom_eval)

void
arb_const_log2(arb_t res, slong prec)
{
    if (prec <= FLINT_BITS * (MP_REAL_CONST_STATIC_N - 1) - 16)
        _arb_const_set_floor(res, MP_REAL_CONST_ID_LOG2, 0, prec);
    else
        _arb_const_round_cached(res, arb_const_log2_hypgeom, prec);
}
