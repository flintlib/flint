/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int TEMPLATE(T, mat_equal) (const TEMPLATE(T, mat_t) mat1,
                            const TEMPLATE(T, mat_t) mat2,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong j;

    if (mat1->r != mat2->r || mat1->c != mat2->c)
    {
        return 0;
    }

    if (mat1->r == 0 || mat1->c == 0)
        return 1;

    for (j = 0; j < mat1->r; j++)
    {
        if (!_TEMPLATE(T, vec_equal)
            (mat1->rows[j], mat2->rows[j], mat1->c, ctx))
        {
            return 0;
        }
    }

    return 1;
}


#endif
