/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_char_dot_slong(ulong a, const slong * n, slong g)
{
    ulong a_shift = a;
    slong sgn = 0;
    slong k;

    for (k = 0; k < g; k++)
    {
        if (a_shift & 1)
        {
            sgn += n[g - 1 - k] & 3;
        }
        a_shift = a_shift >> 1;
    }

    return sgn % 4;
}
