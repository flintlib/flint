/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

ulong
dirichlet_group_num_primitive(const dirichlet_group_t G)
{
    if (G->q % 4 == 2)
        return 0;
    else
    {
        slong k;
        ulong n = 1;

        /* no overflow since result < G->q */

        for (k = (G->neven == 2); k < G->num; k++)
        {
            ulong p = G->P[k].p, e = G->P[k].e;
            if (e == 1)
                n *= p - 2;
            else
                n *= (p * (p - 2) + 1) * n_pow(p, e - 2);
        }

        return n;
    }
}
