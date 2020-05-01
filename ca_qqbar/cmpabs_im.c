/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int
ca_qqbar_cmpabs_im(const ca_qqbar_t x, const ca_qqbar_t y)
{
    int sx, sy;

    if (!arb_overlaps(acb_imagref(CA_QQBAR_ENCLOSURE(x)), acb_imagref(CA_QQBAR_ENCLOSURE(y))))
    {
        return arf_cmpabs(arb_midref(acb_imagref(CA_QQBAR_ENCLOSURE(x))),
                       arb_midref(acb_imagref(CA_QQBAR_ENCLOSURE(y))));
    }

    sx = ca_qqbar_sgn_im(x);
    sy = ca_qqbar_sgn_im(y);

    if (sx == sy)
        return 0;

    if (sy == 0)
        return 1;

    if (sx == 0)
        return -1;

    if (sx > 0 && sy > 0)
        return ca_qqbar_cmp_im(x, y);
    else
        return -ca_qqbar_cmp_im(x, y);
}
