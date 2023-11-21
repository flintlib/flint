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
acb_theta_jet_index(const slong * tup, slong g)
{
    slong ord, res, k;
    slong j;

    /* Get total derivation order */
    ord = acb_theta_jet_total_order(tup, g);
    if (ord == 0 || g == 1)
    {
        return ord;
    }

    /* Count tuples with smaller total order */
    res = acb_theta_jet_nb(ord - 1, g);

    for (j = 0; j < g - 1; j++)
    {
        k = tup[j];
        res += acb_theta_jet_nb(ord - k - 1, g - 1 - j);
        ord -= k;
    }

    return res;
}
