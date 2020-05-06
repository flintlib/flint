/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
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
        flint_printf("qqbar_tan_pi: division by zero\n");
        flint_abort();
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
}
