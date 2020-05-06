/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int
qqbar_cmpabs_re(const qqbar_t x, const qqbar_t y)
{
    int sx, sy;

    if (!arb_overlaps(acb_realref(QQBAR_ENCLOSURE(x)), acb_realref(QQBAR_ENCLOSURE(y))))
    {
        return arf_cmpabs(arb_midref(acb_realref(QQBAR_ENCLOSURE(x))),
                       arb_midref(acb_realref(QQBAR_ENCLOSURE(y))));
    }

    sx = qqbar_sgn_re(x);
    sy = qqbar_sgn_re(y);

    if (sx == sy)
        return 0;

    if (sy == 0)
        return 1;

    if (sx == 0)
        return -1;

    if (sx > 0 && sy > 0)
        return qqbar_cmp_re(x, y);
    else
        return -qqbar_cmp_re(x, y);
}
