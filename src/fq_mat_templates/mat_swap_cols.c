/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_swap_cols)(TEMPLATE(T, mat_t) mat, slong * perm, slong r, slong s, const TEMPLATE(T, ctx_t) ctx)
{
    if (r != s && !TEMPLATE(T, mat_is_empty)(mat, ctx))
    {
        slong t;

        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        for (t = 0; t < mat->r; t++)
            TEMPLATE(T, swap)(TEMPLATE(T, mat_entry)(mat, t, r), TEMPLATE(T, mat_entry)(mat, t, s), ctx);
    }
}

#endif
