/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

slong
fmpzi_canonical_unit_i_pow(const fmpzi_t x)
{
    int s, t;

    s = fmpz_cmp(fmpzi_realref(x), fmpzi_imagref(x));

    if (s == 0)
    {
        if (fmpz_sgn(fmpzi_realref(x)) < 0)
            return 2;
        else
            return 0;
    }
    else
    {
        t = fmpz_cmpabs(fmpzi_realref(x), fmpzi_imagref(x));

        if (s > 0)
            return (t <= 0) ? 1 : 0;
        else
            return (t <= 0) ? 3 : 2;
    }
}
