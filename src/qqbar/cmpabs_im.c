/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

static int
_arb_overlaps_abs(const arb_t x, const arb_t y)
{
    arb_t t, u;

    *t = *x;
    *u = *y;

    if (arf_sgn(arb_midref(t)) < 0)
        ARF_NEG(arb_midref(t));
    if (arf_sgn(arb_midref(u)) < 0)
        ARF_NEG(arb_midref(u));

    return arb_overlaps(t, u);
}

int
qqbar_cmpabs_im(const qqbar_t x, const qqbar_t y)
{
    int sx, sy;

    if (!_arb_overlaps_abs(acb_imagref(QQBAR_ENCLOSURE(x)), acb_imagref(QQBAR_ENCLOSURE(y))))
    {
        return arf_cmpabs(arb_midref(acb_imagref(QQBAR_ENCLOSURE(x))),
                       arb_midref(acb_imagref(QQBAR_ENCLOSURE(y))));
    }

    sx = qqbar_sgn_im(x);
    sy = qqbar_sgn_im(y);

    if (sx == 0 && sy == 0)
        return 0;
    if (sy == 0 && sx != 0)
        return 1;
    if (sx == 0 && sy != 0)
        return -1;
    if (sx > 0 && sy > 0)
        return qqbar_cmp_im(x, y);
    if (sx < 0 && sy < 0)
        return -qqbar_cmp_im(x, y);

    /* opposite signs */
    {
        int res;
        qqbar_t t;
        qqbar_init(t);
        if (sx > 0)
        {
            qqbar_neg(t, y);
            res = qqbar_cmp_im(x, t);
        }
        else
        {
            qqbar_neg(t, x);
            res = qqbar_cmp_im(t, y);
        }
        qqbar_clear(t);
        return res;
    }
}
