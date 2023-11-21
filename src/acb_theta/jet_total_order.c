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
acb_theta_jet_total_order(const slong * tup, slong g)
{
    slong k;
    slong res = 0;

    for (k = 0; k < g; k++)
    {
        res += tup[k];
    }

    return res;
}
