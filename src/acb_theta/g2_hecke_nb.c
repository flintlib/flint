/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_g2_hecke_nb(slong q)
{
    slong p;
    int is_T1;

    if (n_is_prime(q))
    {
        p = q;
        is_T1 = 0;
    }
    else
    {
        p = n_sqrt(q);
        is_T1 = 1;
        if (p * p != q || !n_is_prime(p))
        {
            return 0;
        }
    }

    if (is_T1)
    {
        return p + n_pow(p, 2) + n_pow(p, 3) + n_pow(p, 4);
    }
    else
    {
        return 1 + p + n_pow(p, 2) + n_pow(p, 3);
    }
}
