/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

slong
fmpr_add_si(fmpr_t z, const fmpr_t x, slong y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_si(t, y);
    r = fmpr_add(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

