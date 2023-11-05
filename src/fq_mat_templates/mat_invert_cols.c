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
TEMPLATE(T, mat_invert_cols)(TEMPLATE(T, mat_t) mat, slong * perm, const TEMPLATE(T, ctx_t) ctx)
{
    if (!TEMPLATE(T, mat_is_empty)(mat, ctx))
    {
        slong t, i;
        slong c = mat->c;
        slong k = mat->c/2;

        if (perm != NULL)
            for (i = 0; i < k; i++)
                FLINT_SWAP(slong, perm[i], perm[c - i - 1]);

        for (t = 0; t < mat->r; t++)
            for (i = 0; i < k; i++)
                TEMPLATE(T, swap)(TEMPLATE(T, mat_entry)(mat, t, i), TEMPLATE(T, mat_entry)(mat, t, c - i - 1), ctx);
    }
}

#endif
