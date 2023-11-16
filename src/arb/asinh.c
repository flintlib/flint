/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_asinh(arb_t z, const arb_t x, slong prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(z);
    }
    else
    {
        arb_t t;
        arb_init(t);

        arb_mul(t, x, x, prec + 4);
        arb_sqrt1pm1(t, t, prec + 4);

        if (arf_sgn(arb_midref(x)) >= 0)
        {
            arb_add(t, t, x, prec + 4);
            arb_log1p(z, t, prec);
        }
        else
        {
            arb_sub(t, t, x, prec + 4);
            arb_log1p(z, t, prec);
            arb_neg(z, z);
        }

        arb_clear(t);
    }
}

