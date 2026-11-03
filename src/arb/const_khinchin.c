/*
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mp_real.h"

static void
arb_const_khinchin_eval(arb_t s, slong prec)
{
    _mp_real_const_arb(s, mp_real_const_khinchin, prec);
}

ARB_DEF_CACHED_CONSTANT(arb_const_khinchin, arb_const_khinchin_eval)
