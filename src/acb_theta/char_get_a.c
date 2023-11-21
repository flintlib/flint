/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

ulong
acb_theta_char_get_a(const slong * n, slong g)
{
    slong k;
    ulong a = 0;

    for (k = 0; k < g; k++)
    {
        a *= 2;
        a += ((n[k] % 2) + 2) % 2;
    }

    return a;
}
