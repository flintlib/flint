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
dirichlet_conductor_ui(const dirichlet_group_t G, ulong a)
{
    slong k;
    ulong ap, cond;

    cond = 1;

    for (k = (G->neven == 2); k < G->num; k++)
    {
        ulong p;
        nmod_t pe;
        p = G->P[k].p;
        pe = G->P[k].pe;
        ap = a % pe.n;
        if (ap == 1)
            continue;
        if (p == 2)
        {
            cond = 4;
            if (a % 4 == 3)
                ap = pe.n - ap;
        }
        else
        {
            cond *= p;
            ap = nmod_pow_ui(ap, p - 1, pe);
        }

        while (ap != 1)
        {
            cond *= p;
            ap = nmod_pow_ui(ap, p, pe);
        }

    }

    return cond;
}
