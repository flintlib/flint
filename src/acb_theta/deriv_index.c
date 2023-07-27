/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_deriv_index(const slong* orders, slong g)
{
    slong ord, res, k;
    slong j;

    if (g == 1)
    {
        return 0;
    }

    /* Get total derivation order */
    ord = 0;
    for (j = 0; j < g; j++)
    {
        ord += orders[j];
    }

    res = 0;
    k = orders[0];
    /* First count all tuples with first entry g, ..., k+1 */
    for (j = 0; j < ord - k; j++)
    {
        res += acb_theta_deriv_nb(j, g - 1);
    }
    /* Now it is the same as index as orders[1] in g - 1 variables */
    res += acb_theta_deriv_index(orders + 1, g - 1);

    return res;
}
