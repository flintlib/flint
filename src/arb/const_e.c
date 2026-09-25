/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mp_real.h"

static void
arb_const_e_eval(arb_t s, slong prec)
{
    /* e = sum 1/k! by binary splitting (mp_real/const_e.c) */
    _mp_real_const_arb(s, mp_real_const_e, prec);
}

ARB_DEF_CACHED_CONSTANT(arb_const_e, arb_const_e_eval)
