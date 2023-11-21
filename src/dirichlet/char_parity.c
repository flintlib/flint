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
dirichlet_parity_char(const dirichlet_group_t G, const dirichlet_char_t x)
{
    int k, odd = 0;
    for (k = 0; k < G->num; k++)
    {
        if (k == 1 && G->neven == 2)
            continue;
        if (x->log[k] % 2)
            odd = 1 - odd;
    }
    return odd;
}
