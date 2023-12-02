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
dirichlet_char_lower(dirichlet_char_t y, const dirichlet_group_t H, const dirichlet_char_t x, const dirichlet_group_t G)
{
    slong k, l;

    if (G->q % H->q != 0)
    {
        flint_throw(FLINT_ERROR, "conrey_lower: lower modulus %wu does not divide %wu\n", H->q, G->q);
    }

    for (k = 0, l = 0; k < G->num && l < H->num; k++)
    {
        ulong p = G->P[k].p;
        if (p == H->P[l].p)
        {
            ulong pef = n_pow(G->P[k].p, G->P[k].e - H->P[l].e);
            ulong a = x->log[k];
            if (a % pef)
            {
                    flint_throw(FLINT_ERROR, "conrey_lower: conductor does not divide lower modulus %wu", H->q);
            }
            y->log[l] = a / pef;
            l++;
        }
    }
    _dirichlet_char_exp(y, H);
}
