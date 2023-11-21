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
qqbar_cmpabs(const qqbar_t x, const qqbar_t y)
{
    slong prec;
    acb_t z1, z2;
    arb_t z3, z4;
    int res;

    if (qqbar_sgn_im(x) == 0 && qqbar_sgn_im(y) == 0)
        return qqbar_cmpabs_re(x, y);

    if (qqbar_sgn_re(x) == 0 && qqbar_sgn_re(y) == 0)
        return qqbar_cmpabs_im(x, y);

    {
        acb_init(z1);
        acb_init(z2);
        arb_init(z3);
        arb_init(z4);

        acb_set(z1, QQBAR_ENCLOSURE(x));
        acb_set(z2, QQBAR_ENCLOSURE(y));

        res = 0;
        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _qqbar_enclosure_raw(z1, QQBAR_POLY(x), z1, prec);
            _qqbar_enclosure_raw(z2, QQBAR_POLY(y), z2, prec);

            acb_abs(z3, z1, prec);
            acb_abs(z4, z2, prec);

            if (!arb_overlaps(z3, z4))
            {
                res = arf_cmpabs(arb_midref(z3), arb_midref(z4));
                break;
            }

            /* Force an exact computation (may be slow) */
            if (prec >= 4 * QQBAR_DEFAULT_PREC)
            {
                qqbar_t t, u;
                qqbar_init(t);
                qqbar_init(u);
                qqbar_abs2(t, x);
                qqbar_abs2(u, y);
                res = qqbar_cmp_re(t, u);
                qqbar_clear(t);
                qqbar_clear(u);
                break;
            }
        }

        acb_clear(z1);
        acb_clear(z2);
        arb_clear(z3);
        arb_clear(z4);
    }

    return res;
}

