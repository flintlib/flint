/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
acb_theta_is_balanced(slong * j0, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_t test;
    slong j;
    int r = 1;

    arb_init(test);

    for (j = 0; j < g - 1; j++)
    {
        arb_mul_si(test, acb_imagref(arb_mat_entry(tau, j, j)),
                   ACB_THETA_BALANCE_THRESHOLD, prec);
        if (arb_lt(test, acb_imagref(arb_mat_entry(tau, j + 1, j + 1))))
        {
            r = 0;
            *j0 = j;
            break;
        }
    }

    arb_clear(test);
    return r;
}
