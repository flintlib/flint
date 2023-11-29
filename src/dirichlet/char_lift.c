/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

void
dirichlet_char_lift(dirichlet_char_t y, const dirichlet_group_t G, const dirichlet_char_t x, const dirichlet_group_t H)
{
    slong k, l;

    if (G->q % H->q != 0)
    {
        flint_throw(FLINT_ERROR, "conrey_lift: lower modulus %wu does not divide %wu\n", H->q, G->q);
    }

    for (k = 0, l = 0; k < G->num && l < H->num; k++)
    {
        if (G->P[k].p == H->P[l].p)
        {
            slong e = n_pow(G->P[k].p, G->P[k].e - H->P[l].e);
            y->log[k] = x->log[l] * e;
            l++;
        }
    }

    _dirichlet_char_exp(y, G);
}
