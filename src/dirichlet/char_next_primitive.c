/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

int
dirichlet_char_next_primitive(dirichlet_char_t x, const dirichlet_group_t G)
{
    /* update index avoiding multiples of p except for first component
       if 8|q */
    slong k;
    /*
    if (G->neven == 2)
    {
        x->n = nmod_mul(x->n, G->generators[0], G->mod);
        x->log[0]++;
        if (x->log[0] == 1)
            return 0;
        x->log[0] = 0;
        k = 1;
    }
    */
    for (k = G->num - 1; k >= 0; k--)
    {
#if 1
        x->n = nmod_mul(x->n, G->generators[k], G->mod);
        x->log[k]++;
        if (x->log[k] % G->P[k].p == 0)
        {
            x->n = nmod_mul(x->n, G->generators[k], G->mod);
            x->log[k]++;
        }
        if (x->log[k] < G->P[k].phi.n)
            break;
        if (x->log[k] == G->P[k].phi.n)
            x->n = nmod_mul(x->n, G->generators[k], G->mod);
        x->log[k] =  1;
#else
        do {
            x->n = nmod_mul(x->n, G->generators[k], G->mod);
            x->log[k]++;
        } while (x->log[k] % G->P[k].p == 0);
        if (x->log[k] < G->P[k].phi)
            break;
        if (x->log[k] == G->P[k].phi)
            x->n = nmod_mul(x->n, G->generators[k], G->mod);
        x->log[k] =  1;
#endif
    }
    /* return last index modified */
    return k;
}
