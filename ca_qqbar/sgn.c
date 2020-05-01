/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_sgn(ca_qqbar_t res, const ca_qqbar_t x)
{
    int re, im;

    re = ca_qqbar_sgn_re(x);
    im = ca_qqbar_sgn_im(x);

    if (im == 0)
    {
        ca_qqbar_set_si(res, re);
    }
    else if (re == 0)
    {
        ca_qqbar_i(res);
        if (im < 0)
            ca_qqbar_neg(res, res);
    }
    else
    {
        ca_qqbar_t t;
        ca_qqbar_init(t);
        ca_qqbar_abs(t, x);
        ca_qqbar_div(res, x, t);
        ca_qqbar_clear(t);
    }
}
