/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_randtest) (TEMPLATE(T, mat_t) mat, flint_rand_t state,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong r, c, i, j;

    r = mat->r;
    c = mat->c;

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            TEMPLATE(T, randtest) (mat->rows[i] + j, state, ctx);
}


#endif
