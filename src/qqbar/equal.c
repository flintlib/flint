/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qqbar.h"

int
qqbar_equal(const qqbar_t x, const qqbar_t y)
{
    slong prec;
    acb_t z1, z2, z3;
    int res;

    if (x == y)
        return 1;

    if (!fmpz_poly_equal(QQBAR_POLY(x), QQBAR_POLY(y)))
        return 0;

    if (qqbar_degree(x) == 1)
        return 1;

    if (!acb_overlaps(QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(y)))
        return 0;

    if (acb_contains(QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(y)))
        return 1;

    if (acb_contains(QQBAR_ENCLOSURE(y), QQBAR_ENCLOSURE(x)))
        return 1;

    acb_init(z1);
    acb_init(z2);
    acb_init(z3);

    acb_set(z1, QQBAR_ENCLOSURE(x));
    acb_set(z2, QQBAR_ENCLOSURE(y));

    res = 0;
    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        _qqbar_enclosure_raw(z1, QQBAR_POLY(x), z1, prec);
        _qqbar_enclosure_raw(z2, QQBAR_POLY(y), z2, prec);

        if (!acb_overlaps(z1, z2))
        {
            res = 0;
            break;
        }

        acb_union(z3, z1, z2, prec);

        if (_qqbar_validate_uniqueness(z3, QQBAR_POLY(x), z3, 2 * prec))
        {
            res = 1;
            break;
        }
    }

    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);

    return res;
}

