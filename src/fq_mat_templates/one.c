/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_one)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, n;

    TEMPLATE(T, mat_zero)(mat, ctx);
    n = FLINT_MIN(mat->r, mat->c);

    for (i = 0; i < n; i++)
        TEMPLATE(T, one)(TEMPLATE(T, mat_entry)(mat, i, i), ctx);
}

#endif

