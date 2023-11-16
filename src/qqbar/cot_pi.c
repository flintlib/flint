/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "qqbar.h"

int
qqbar_cot_pi(qqbar_t res, slong p, ulong q)
{
    slong g;

    g = n_gcd(FLINT_ABS(p), q);

    if (g != 1)
    {
        p /= g;
        q /= g;
    }

    if (q == 1)
    {
        return 0;
    }
    else if (q == 2)
    {
        qqbar_zero(res);
    }
    else
    {
        qqbar_tan_pi(res, p, q);
        qqbar_inv(res, res);
    }

    return 1;
}
