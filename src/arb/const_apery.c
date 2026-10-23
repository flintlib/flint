/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mp_real.h"

static void
arb_const_apery_eval(arb_t s, slong prec)
{
    /* Zuniga's 2023-vi series (fixed/const_zeta3.c) */
    _mp_real_const_arb(s, mp_real_const_zeta3, prec);
}

ARB_DEF_CACHED_CONSTANT(arb_const_apery, arb_const_apery_eval)
