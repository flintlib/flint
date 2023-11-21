/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_siegel_randtest_vec(acb_ptr z, flint_rand_t state, slong g, slong prec)
{
    slong mag_bits = n_randint(state, 4);
    slong k;

    for (k = 0; k < g; k++)
    {
        switch (n_randint(state, 10))
        {
        case 0:
            acb_randtest_param(&z[k], state, prec, mag_bits);
            break;
        case 1:
            acb_randtest(&z[k], state, prec, mag_bits);
            break;
        case 2:
            acb_randtest_precise(&z[k], state, prec, mag_bits);
            break;
        case 3:
            acb_randtest(&z[k], state, prec, 20);
            break;
        default:
            acb_urandom(&z[k], state, prec);
        }
    }
}
