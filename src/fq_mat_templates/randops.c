/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_randops) (TEMPLATE(T, mat_t) mat,
                          flint_rand_t state, slong count, const TEMPLATE(T, ctx_t) ctx)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;

    if (mat->r == 0 || mat->c == 0)
        return;

    for (c = 0; c < count; c++)
    {
        if (n_randint(state, 2))
        {
            if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
                continue;

            if (n_randint(state, 2))
                for (k = 0; k < n; k++)
                    TEMPLATE(T, add) (TEMPLATE(T, mat_entry)(mat, j, k),
                                      TEMPLATE(T, mat_entry)(mat, j, k), TEMPLATE(T, mat_entry)(mat, i, k), ctx);
            else
                for (k = 0; k < n; k++)
                    TEMPLATE(T, sub) (TEMPLATE(T, mat_entry)(mat, j, k),
                                      TEMPLATE(T, mat_entry)(mat, j, k), TEMPLATE(T, mat_entry)(mat, i, k), ctx);
        }
        else
        {
            if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                continue;
            if (n_randint(state, 2))
                for (k = 0; k < m; k++)
                    TEMPLATE(T, add) (TEMPLATE(T, mat_entry)(mat, k, j),
                                      TEMPLATE(T, mat_entry)(mat, k, j), TEMPLATE(T, mat_entry)(mat, k, i), ctx);
            else
                for (k = 0; k < m; k++)
                    TEMPLATE(T, sub) (TEMPLATE(T, mat_entry)(mat, k, j),
                                      TEMPLATE(T, mat_entry)(mat, k, j), TEMPLATE(T, mat_entry)(mat, k, i), ctx);
        }
    }
}


#endif
