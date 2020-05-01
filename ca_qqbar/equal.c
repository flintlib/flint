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
ca_qqbar_equal(const ca_qqbar_t x, const ca_qqbar_t y)
{
    slong prec;
    acb_t z1, z2, z3;
    int res;

    if (x == y)
        return 1;

    if (!fmpz_poly_equal(CA_QQBAR_POLY(x), CA_QQBAR_POLY(y)))
        return 0;

    if (ca_qqbar_degree(x) == 1)
        return 1;

    if (!acb_overlaps(CA_QQBAR_ENCLOSURE(x), CA_QQBAR_ENCLOSURE(y)))
        return 0;

    if (acb_contains(CA_QQBAR_ENCLOSURE(x), CA_QQBAR_ENCLOSURE(y)))
        return 1;

    if (acb_contains(CA_QQBAR_ENCLOSURE(y), CA_QQBAR_ENCLOSURE(x)))
        return 1;

    acb_init(z1);
    acb_init(z2);
    acb_init(z3);

    acb_set(z1, CA_QQBAR_ENCLOSURE(x));
    acb_set(z2, CA_QQBAR_ENCLOSURE(y));

    res = 0;
    for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        _ca_qqbar_enclosure_raw(z1, CA_QQBAR_POLY(x), z1, prec);
        _ca_qqbar_enclosure_raw(z2, CA_QQBAR_POLY(y), z2, prec);

        if (!acb_overlaps(z1, z2))
        {
            res = 0;
            break;
        }

        acb_union(z3, z1, z2, prec);

        if (_ca_qqbar_validate_enclosure(z3, CA_QQBAR_POLY(x), z3, 2 * prec))
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

