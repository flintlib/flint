/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
_arb_const_log_sqrt2pi(arb_t t, slong prec)
{
    arb_const_pi(t, prec + 2);
    arb_mul_2exp_si(t, t, 1);
    arb_log(t, t, prec);
    arb_mul_2exp_si(t, t, -1);
}

ARB_DEF_CACHED_CONSTANT(arb_const_log_sqrt2pi, _arb_const_log_sqrt2pi)
