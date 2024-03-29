/*
    Copyright (C) 2010,2011 Fredrik Johansson
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
TEMPLATE(T, mat_solve_triu_classical) (TEMPLATE(T, mat_t) X,
                                       const TEMPLATE(T, mat_t) U,
                                       const TEMPLATE(T, mat_t) B,
                                       int unit, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, n, m;
    TEMPLATE(T, struct) * inv, *tmp;

    n = U->r;
    m = B->c;

    if (!unit)
    {
        inv = _TEMPLATE(T, vec_init) (n, ctx);
        for (i = 0; i < n; i++)
            TEMPLATE(T, inv) (inv + i, TEMPLATE(T, mat_entry) (U, i, i), ctx);
    }
    else
        inv = NULL;

    tmp = _TEMPLATE(T, vec_init) (n, ctx);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            TEMPLATE(T, set) (tmp + j, TEMPLATE(T, mat_entry) (X, j, i), ctx);

        for (j = n - 1; j >= 0; j--)
        {
            TEMPLATE(T, t) s;
            TEMPLATE(T, init) (s, ctx);
            _TEMPLATE(T, vec_dot) (s, U->rows[j] + j + 1, tmp + j + 1,
                                   n - j - 1, ctx);
            TEMPLATE(T, sub) (s, TEMPLATE(T, mat_entry) (B, j, i), s, ctx);
            if (!unit)
                TEMPLATE(T, mul) (s, s, inv + j, ctx);
            TEMPLATE(T, set) (tmp + j, s, ctx);
            TEMPLATE(T, clear) (s, ctx);
        }

        for (j = 0; j < n; j++)
            TEMPLATE(T, mat_entry_set) (X, j, i, tmp + j, ctx);
    }

    _TEMPLATE(T, vec_clear) (tmp, n, ctx);
    if (!unit)
        _TEMPLATE(T, vec_clear) (inv, n, ctx);
}


#endif
