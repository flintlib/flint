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
dirichlet_conductor_char(const dirichlet_group_t G, const dirichlet_char_t x)
{
    int k, f;
    ulong cond = 1;

    if (G->neven >= 1 && x->log[0] == 1)
        cond = 4;
    if (G->neven == 2 && x->log[1])
    {
        ulong l2 = x->log[1];
        f = n_remove(&l2, 2);
        cond = UWORD(1) << (G->P[1].e - f);
    }

    for (k = G->neven; k < G->num; k++)
    {
        if (x->log[k])
        {
            ulong p, lp;
            p = G->P[k].p;
            lp = x->log[k];
            f = n_remove(&lp, p);
            if (f)
                cond *= n_pow(p, G->P[k].e - f);
            else
                cond *= G->P[k].pe.n;
        }
    }

    return cond;
}
