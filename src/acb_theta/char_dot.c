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
acb_theta_char_dot(ulong a, ulong b, slong g)
{
    int sgn = 0;
    slong k;
    ulong and = a & b;

    for (k = 0; k < g; k++)
    {
        if (and & 1)
        {
            sgn++;
        }
        and = and >> 1;
    }
    return sgn % 4;
}
