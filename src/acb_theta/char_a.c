/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

ulong
acb_theta_char_a(slong* coords, slong g)
{
    ulong a = 0;
    slong k;

    for (k = 0; k < g; k++)
    {
        a = a << 1;
        a += ((4 + coords[k] % 4) % 4) / 2;
    }

    return a;
}
