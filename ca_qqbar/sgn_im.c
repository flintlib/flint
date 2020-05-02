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
ca_qqbar_sgn_im(const ca_qqbar_t x)
{
    if (ca_qqbar_degree(x) == 1)
    {
        return 0;
    }
    else if (arb_is_zero(acb_imagref(CA_QQBAR_ENCLOSURE(x))))
    {
        return 0;
    }
    else if (!arb_contains_zero(acb_imagref(CA_QQBAR_ENCLOSURE(x))))
    {
        return arf_sgn(arb_midref(acb_imagref(CA_QQBAR_ENCLOSURE(x))));
    }
    else
    {
        slong prec;
        int res;
        acb_t t, u;
        acb_init(t);
        acb_init(u);

        acb_set(t, CA_QQBAR_ENCLOSURE(x));
        res = 0;

        for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _ca_qqbar_enclosure_raw(t, CA_QQBAR_POLY(x), t, prec);

            if (!arb_contains_zero(acb_imagref(t)))
            {
                res = arf_sgn(arb_midref(acb_imagref(t)));
                break;
            }

            acb_conj(u, t);
            acb_union(u, u, t, prec);

            if (_ca_qqbar_validate_uniqueness(u, CA_QQBAR_POLY(x), u, 2 * prec))
                break;
        }

        acb_clear(t);
        acb_clear(u);

        return res;
    }
}
