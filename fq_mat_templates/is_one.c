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

int
TEMPLATE(T, mat_is_one) (const TEMPLATE(T, mat_t) mat,
                          const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
            if (i == j ? !TEMPLATE(T, is_one)(TEMPLATE(T, mat_entry)(mat, i, j), ctx) :
                         !TEMPLATE(T, is_zero)(TEMPLATE(T, mat_entry)(mat, i, j), ctx))
                return 0;
    }

    return 1;
}

#endif
