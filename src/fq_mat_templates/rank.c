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

slong
TEMPLATE(T, mat_rank) (const TEMPLATE(T, mat_t) A,
                       const TEMPLATE(T, ctx_t) ctx)
{
    slong m, n, rank;
    slong *perm;
    TEMPLATE(T, mat_t) tmp;

    m = A->r;
    n = A->c;

    if (m == 0 || n == 0)
        return 0;

    TEMPLATE(T, mat_init_set) (tmp, A, ctx);
    perm = flint_malloc(sizeof(slong) * m);

    rank = TEMPLATE(T, mat_lu) (perm, tmp, 0, ctx);

    flint_free(perm);
    TEMPLATE(T, mat_clear) (tmp, ctx);
    return rank;
}


#endif
