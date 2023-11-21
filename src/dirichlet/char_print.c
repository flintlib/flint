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
dirichlet_char_print(const dirichlet_group_t G, const dirichlet_char_t x)
{
    slong k;
    if (G->num)
        flint_printf("[%wu", x->log[0]);
    else
        flint_printf("[");
    for (k = 1; k < G->num; k++)
        flint_printf(", %wu", x->log[k]);
    flint_printf("]");
}
