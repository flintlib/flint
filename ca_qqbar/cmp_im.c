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
ca_qqbar_cmp_im(const ca_qqbar_t x, const ca_qqbar_t y)
{
    slong prec;
    acb_t z1, z2;
    int res;

    if (!arb_overlaps(acb_imagref(CA_QQBAR_ENCLOSURE(x)), acb_imagref(CA_QQBAR_ENCLOSURE(y))))
    {
        return arf_cmp(arb_midref(acb_imagref(CA_QQBAR_ENCLOSURE(x))),
                       arb_midref(acb_imagref(CA_QQBAR_ENCLOSURE(y))));
    }

    if (ca_qqbar_imag_sgn(y) == 0)
        return ca_qqbar_imag_sgn(x);

    if (ca_qqbar_imag_sgn(x) == 0)
        return -ca_qqbar_imag_sgn(y);

    if (ca_qqbar_equal(x, y))
        return 0;

    {
        acb_init(z1);
        acb_init(z2);

        acb_set(z1, CA_QQBAR_ENCLOSURE(x));
        acb_set(z2, CA_QQBAR_ENCLOSURE(y));

        res = 0;
        for (prec = CA_QQBAR_DEFAULT_PREC; ; prec *= 2)
        {
            _ca_qqbar_enclosure_raw(z1, CA_QQBAR_POLY(x), z1, prec);
            _ca_qqbar_enclosure_raw(z2, CA_QQBAR_POLY(y), z2, prec);

            if (!arb_overlaps(acb_imagref(z1), acb_imagref(z2)))
            {
                res = arf_cmp(arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)));
                break;
            }

            /* Force an exact computation (may be slow) */
            if (prec >= 4 * CA_QQBAR_DEFAULT_PREC)
            {
                ca_qqbar_t t;
                ca_qqbar_init(t);
                ca_qqbar_sub(t, x, y);
                res = ca_qqbar_imag_sgn(t);
                ca_qqbar_clear(t);
                break;
            }
        }

        acb_clear(z1);
        acb_clear(z2);
    }

    return res;
}

