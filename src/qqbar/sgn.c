/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_sgn(qqbar_t res, const qqbar_t x)
{
    int re, im;

    re = qqbar_sgn_re(x);
    im = qqbar_sgn_im(x);

    if (im == 0)
    {
        qqbar_set_si(res, re);
    }
    else if (re == 0)
    {
        qqbar_i(res);
        if (im < 0)
            qqbar_neg(res, res);
    }
    else
    {
        qqbar_t t;
        qqbar_init(t);
        qqbar_abs(t, x);
        qqbar_div(res, x, t);
        qqbar_clear(t);
    }
}
