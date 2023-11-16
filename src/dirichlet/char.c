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
dirichlet_char_init(dirichlet_char_t x, const dirichlet_group_t G) {
    x->log = flint_malloc(G->num * sizeof(ulong));
    dirichlet_char_one(x, G);
}

void
dirichlet_char_clear(dirichlet_char_t x) {
    flint_free(x->log);
}
