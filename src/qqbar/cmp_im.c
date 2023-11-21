/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int
qqbar_cmp_im(const qqbar_t x, const qqbar_t y)
{
    slong prec;
    acb_t z1, z2;
    int res;

    if (!arb_overlaps(acb_imagref(QQBAR_ENCLOSURE(x)), acb_imagref(QQBAR_ENCLOSURE(y))))
    {
        return arf_cmp(arb_midref(acb_imagref(QQBAR_ENCLOSURE(x))),
                       arb_midref(acb_imagref(QQBAR_ENCLOSURE(y))));
    }

    if (qqbar_sgn_im(y) == 0)
        return qqbar_sgn_im(x);

    if (qqbar_sgn_im(x) == 0)
        return -qqbar_sgn_im(y);

    if (qqbar_equal(x, y))
        return 0;

    {
        qqbar_t t;
        qqbar_init(t);
        qqbar_neg(t, y);
        qqbar_conj(t, t);
        res = qqbar_equal(x, t);
        qqbar_clear(t);
        if (res == 1)
            return 0;
    }

    acb_init(z1);
    acb_init(z2);

    acb_set(z1, QQBAR_ENCLOSURE(x));
    acb_set(z2, QQBAR_ENCLOSURE(y));

    res = 0;
    for (prec = QQBAR_DEFAULT_PREC; ; prec *= 2)
    {
        _qqbar_enclosure_raw(z1, QQBAR_POLY(x), z1, prec);
        _qqbar_enclosure_raw(z2, QQBAR_POLY(y), z2, prec);

        if (!arb_overlaps(acb_imagref(z1), acb_imagref(z2)))
        {
            res = arf_cmp(arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)));
            break;
        }

        /* Force an exact computation (may be slow) */
        if (prec >= 4 * QQBAR_DEFAULT_PREC)
        {
            qqbar_t t;
            qqbar_init(t);
            qqbar_sub(t, x, y);
            res = qqbar_sgn_im(t);
            qqbar_clear(t);
            break;
        }
    }

    acb_clear(z1);
    acb_clear(z2);

    return res;
}

