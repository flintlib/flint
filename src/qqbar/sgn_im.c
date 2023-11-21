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
qqbar_sgn_im(const qqbar_t x)
{
    if (qqbar_degree(x) == 1)
    {
        return 0;
    }
    else if (arb_is_zero(acb_imagref(QQBAR_ENCLOSURE(x))))
    {
        return 0;
    }
    else if (!arb_contains_zero(acb_imagref(QQBAR_ENCLOSURE(x))))
    {
        return arf_sgn(arb_midref(acb_imagref(QQBAR_ENCLOSURE(x))));
    }
    else
    {
        slong prec;
        int res;
        acb_t t, u;
        acb_init(t);
        acb_init(u);

        acb_set(t, QQBAR_ENCLOSURE(x));
        res = 0;

        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _qqbar_enclosure_raw(t, QQBAR_POLY(x), t, prec);

            if (!arb_contains_zero(acb_imagref(t)) || arb_is_zero(acb_imagref(t)))
            {
                res = arf_sgn(arb_midref(acb_imagref(t)));
                break;
            }

#if 0
            acb_conj(u, t);
            acb_union(u, u, t, prec);

            if (_qqbar_validate_uniqueness(u, QQBAR_POLY(x), u, 2 * prec))
                break;
#else
            acb_set(u, t);
            arb_zero(acb_imagref(u));

            if (_qqbar_validate_existence_uniqueness(u, QQBAR_POLY(x), u, 2 * prec))
                break;
#endif
        }

        acb_clear(t);
        acb_clear(u);

        return res;
    }
}
