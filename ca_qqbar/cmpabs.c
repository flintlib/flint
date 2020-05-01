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
ca_qqbar_cmpabs(const ca_qqbar_t x, const ca_qqbar_t y)
{
    slong prec;
    acb_t z1, z2;
    arb_t z3, z4;
    int res;

    if (ca_qqbar_imag_sgn(x) == 0 && ca_qqbar_imag_sgn(y) == 0)
        return ca_qqbar_cmpabs_re(x, y);

    if (ca_qqbar_real_sgn(x) == 0 && ca_qqbar_real_sgn(y) == 0)
        return ca_qqbar_cmpabs_im(x, y);

    {
        acb_init(z1);
        acb_init(z2);
        arb_init(z3);
        arb_init(z4);

        acb_set(z1, CA_QQBAR_ENCLOSURE(x));
        acb_set(z2, CA_QQBAR_ENCLOSURE(y));

        res = 0;
        for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _ca_qqbar_enclosure_raw(z1, CA_QQBAR_POLY(x), z1, prec);
            _ca_qqbar_enclosure_raw(z2, CA_QQBAR_POLY(y), z2, prec);

            acb_abs(z3, z1, prec);
            acb_abs(z4, z2, prec);

            if (!arb_overlaps(z3, z4))
            {
                res = arf_cmpabs(arb_midref(z3), arb_midref(z4));
                break;
            }

            /* Force an exact computation (may be slow) */
            if (prec >= 4 * CA_QQBAR_DEFAULT_PREC)
            {
                ca_qqbar_t t, u;
                ca_qqbar_init(t);
                ca_qqbar_init(u);
                ca_qqbar_abs2(t, x);
                ca_qqbar_abs2(u, y);
                res = ca_qqbar_cmp_re(t, u);
                ca_qqbar_clear(t);
                ca_qqbar_clear(u);
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

