/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

int
fmpr_cmpabs(const fmpr_t x, const fmpr_t y)
{
    int res, xsign, ysign;
    fmpr_t t;

    if (fmpr_equal(x, y))
        return 0;

    if (fmpr_is_special(x) || fmpr_is_special(y))
    {
        if (fmpr_is_nan(x) || fmpr_is_nan(y))
            return 0;
        if (fmpr_is_zero(x)) return -1;
        if (fmpr_is_zero(y)) return 1;
        if (fmpr_is_inf(x)) return fmpr_is_inf(y) ? 0 : 1;
        if (fmpr_is_inf(y)) return -1;
        return -1;
    }

    /* Reduces to integer comparison if bottom exponents are the same */
    if (fmpz_equal(fmpr_expref(x), fmpr_expref(y)))
    {
        res = fmpz_cmpabs(fmpr_manref(x), fmpr_manref(y));
        if (res != 0)
            res = (res < 0) ? -1 : 1;
    }
    else
    {
        /* TODO: compare position of top exponents to avoid subtraction */
        xsign = fmpr_sgn(x);
        ysign = fmpr_sgn(y);

        fmpr_init(t);
        if (xsign == ysign)
            fmpr_sub(t, x, y, 2, FMPR_RND_DOWN);
        else
            fmpr_add(t, x, y, 2, FMPR_RND_DOWN);
        res = fmpr_sgn(t) * xsign;
        fmpr_clear(t);
    }

    return res;
}
