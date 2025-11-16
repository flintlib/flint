/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2009 Thomas Boothby
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"

static int n_is_square_and_get_sqrt(ulong * s, ulong x)
{
    ulong sq = sqrt((double) x) + 0.5;

    *s = sq;
    return x == sq * sq;
}

#define FLINT_ONE_LINE_MULTIPLIER 480

ulong n_factor_one_line(ulong n, ulong iters)
{
    ulong orig_n = n, in, square, sqrti, mod, factor, factoring = iters, iin;
    n *= FLINT_ONE_LINE_MULTIPLIER;

    iin = 0;
    in = n;
    while (factoring && (iin < in))
    {
        sqrti = n_sqrt(in);
        sqrti++;
        square = sqrti*sqrti;
        mod = square - in;
        if (n_is_square_and_get_sqrt(&factor, mod))
        {
            sqrti -= factor;
            factor = n_gcd(orig_n, sqrti);
            if (factor != UWORD(1))
            {
                return factor;
            }
        }
        factoring--;
        iin = in;
        in += n;
    }

    return 0;
}
