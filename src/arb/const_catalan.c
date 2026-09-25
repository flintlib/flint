/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "fixed.h"

static void
arb_const_catalan_eval(arb_t s, slong prec)
{
    /* Pilehrood's short series (fixed/const_catalan.c) */
    _fixed_constant_arb(s, fball_const_catalan, prec);
}

ARB_DEF_CACHED_CONSTANT(arb_const_catalan, arb_const_catalan_eval)
