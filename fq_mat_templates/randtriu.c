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

void
TEMPLATE(T, mat_randtriu) (TEMPLATE(T, mat_t) mat, flint_rand_t state,
                           int unit, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * e;
    slong i, j;

    for (i = 0; i < mat->r; i++)
    {
        for (j = 0; j < mat->c; j++)
        {
            e = TEMPLATE(T, mat_entry) (mat, i, j);
            if (j > i)
            {
                TEMPLATE(T, randtest) (e, state, ctx);
            }
            else if (i == j)
            {
                TEMPLATE(T, randtest) (e, state, ctx);
                if (unit || TEMPLATE(T, is_zero)(e, ctx))
                    TEMPLATE(T, one) (e, ctx);
            }
            else
            {
                TEMPLATE(T, zero) (e, ctx);
            }
        }
    }
}


#endif
