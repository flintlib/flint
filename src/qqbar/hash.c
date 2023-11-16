/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "calcium.h"
#include "qqbar.h"

ulong qqbar_hash(const qqbar_t x)
{
    ulong s;
    slong i, len;
    const fmpz * c;

    len = QQBAR_POLY(x)->length;
    c = QQBAR_COEFFS(x);

    s = 1234567;

    for (i = 0; i < len; i++)
        s = calcium_fmpz_hash(c + i) * 1000003 + s;

    /* todo: add in some bits describing the enclosure, e.g.
       a floor value or just the signs */

    return s;
}
