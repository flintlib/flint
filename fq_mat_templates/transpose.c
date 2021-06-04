/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_transpose) (TEMPLATE(T, mat_t) B,
                            const TEMPLATE(T, mat_t) A,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    TEMPLATE(T, t) tmp;

    if (A == B)
    {
        TEMPLATE(T, init) (tmp, ctx);
        for (i = 0; i < A->r; ++i)
        {
            for (j = i+1; j < A->c; ++j)
            {
                TEMPLATE(T, set) (tmp, &A->rows[i][j], ctx);
                TEMPLATE(T, set) (&B->rows[i][j], &A->rows[j][i], ctx);
                TEMPLATE(T, set) (&B->rows[j][i], tmp, ctx);
            }
        }
        TEMPLATE(T, clear) (tmp, ctx);
    }
    else
    {
        for (i = 0; i < A->r; ++i)
        {
            for (j = 0; j < A->c; ++j)
            {
                TEMPLATE(T, set) (&B->rows[j][i], &A->rows[i][j], ctx);
            }
        }
    }
    
}


#endif
