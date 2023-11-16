/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_roots_clear(acb_dirichlet_roots_t t)
{
    slong k;

    for (k = 0; k < t->depth; k++)
        _acb_vec_clear(t->Z[k], t->size + 1);

    flint_free(t->Z);
    acb_clear(t->z);
}

