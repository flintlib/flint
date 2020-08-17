/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
TEMPLATE(T, mat_is_zero) (const TEMPLATE(T, mat_t) mat,
                          const TEMPLATE(T, ctx_t) ctx)
{
    slong j;

    if (mat->r == 0 || mat->c == 0)
        return 1;

    for (j = 0; j < mat->r; j++)
    {
        if (!_TEMPLATE(T, vec_is_zero) (mat->rows[j], mat->c, ctx))
            return 0;
    }

    return 1;
}


#endif
