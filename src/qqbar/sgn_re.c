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
qqbar_sgn_re(const qqbar_t x)
{
    if (qqbar_degree(x) == 1)
    {
        return -fmpz_sgn(QQBAR_COEFFS(x));
    }
    else if (arb_is_zero(acb_realref(QQBAR_ENCLOSURE(x))))
    {
        return 0;
    }
    else if (!arb_contains_zero(acb_realref(QQBAR_ENCLOSURE(x))))
    {
        return arf_sgn(arb_midref(acb_realref(QQBAR_ENCLOSURE(x))));
    }
    else
    {
        slong d, i;
        slong prec;
        int res, maybe_zero;
        acb_t t, u;
        acb_init(t);
        acb_init(u);

        d = qqbar_degree(x);

        maybe_zero = 1;
        for (i = 1; i < d && maybe_zero; i += 2)
            if (!fmpz_is_zero(QQBAR_COEFFS(x) + i))
                maybe_zero = 0;

        acb_set(t, QQBAR_ENCLOSURE(x));
        res = 0;

        for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
        {
            _qqbar_enclosure_raw(t, QQBAR_POLY(x), t, prec);

            if (!arb_contains_zero(acb_realref(t)) || arb_is_zero(acb_realref(t)))
            {
                res = arf_sgn(arb_midref(acb_realref(t)));
                break;
            }

            if (maybe_zero)
            {
                acb_set(u, t);
                arb_zero(acb_realref(u));
                if (_qqbar_validate_existence_uniqueness(u, QQBAR_POLY(x), u, prec * 2))
                {
                    res = 0;
                    break;
                }
            }
        }

        acb_clear(t);
        acb_clear(u);

        return res;

    }
}
