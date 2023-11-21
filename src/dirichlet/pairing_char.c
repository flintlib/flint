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
dirichlet_pairing_char(const dirichlet_group_t G, const dirichlet_char_t chi, const dirichlet_char_t x)
{
    ulong v = 0, k;
    nmod_t order;

    nmod_init(&order, G->expo);

    for (k = 0; k < G->num; k++)
        v = nmod_add(v, G->PHI[k] * nmod_mul(chi->log[k], x->log[k], G->P[k].phi), order);

    return v;
}
