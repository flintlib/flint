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

void
TEMPLATE(T, mat_randrank) (TEMPLATE(T, mat_t) mat, flint_rand_t state,
                           slong rank, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, struct) * diag;

    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_throw(FLINT_ERROR, "(%s): Impossible rank.\n", __func__);
    }

    diag = _TEMPLATE(T, vec_init) (rank, ctx);
    for (i = 0; i < rank; i++)
        TEMPLATE(T, randtest_not_zero) (diag + i, state, ctx);

    TEMPLATE(T, mat_randpermdiag) (mat, state, diag, rank, ctx);

    _TEMPLATE(T, vec_clear) (diag, rank, ctx);
}


#endif
