/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_submul(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    if (arb_is_zero(acb_imagref(y)))
    {
        arb_submul(acb_imagref(z), acb_imagref(x), acb_realref(y), prec);
        arb_submul(acb_realref(z), acb_realref(x), acb_realref(y), prec);
    }
    else if (arb_is_zero(acb_imagref(x)))
    {
        arb_submul(acb_imagref(z), acb_imagref(y), acb_realref(x), prec);
        arb_submul(acb_realref(z), acb_realref(y), acb_realref(x), prec);
    }
    else
    {
        acb_t t;
        acb_init(t);
        acb_mul(t, x, y, prec);
        acb_sub(z, z, t, prec);
        acb_clear(t);
    }
}

