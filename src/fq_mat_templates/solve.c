/*
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "templates.h"

int
TEMPLATE(T, mat_solve)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) A,
                           const TEMPLATE(T, mat_t) B, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, rank, *perm;
    TEMPLATE(T, mat_t) LU;
    int result;

    if (A->r == 0 || B->c == 0)
        return 1;

    TEMPLATE(T, mat_init_set)(LU, A, ctx);
    perm = flint_malloc(sizeof(slong) * A->r);
    for (i = 0; i < A->r; i++)
        perm[i] = i;

    rank = TEMPLATE(T, mat_lu)(perm, LU, 1, ctx);

    if (rank == A->r)
    {
        TEMPLATE(T, mat_t) PB;
        TEMPLATE(T, mat_window_init)(PB, B, 0, 0, B->r, B->c, ctx);
        for (i = 0; i < A->r; i++)
            PB->rows[i] = B->rows[perm[i]];

        TEMPLATE(T, mat_solve_tril)(X, LU, PB, 1, ctx);
        TEMPLATE(T, mat_solve_triu)(X, LU, X, 0, ctx);

        TEMPLATE(T, mat_window_clear)(PB, ctx);
        result = 1;
    }
    else
    {
        result = 0;
    }

    TEMPLATE(T, mat_clear)(LU, ctx);
    flint_free(perm);

    return result;
}

#endif
