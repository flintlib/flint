/*
    Copyright (C) 2010,2011 Fredrik Johansson
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
TEMPLATE(T, mat_solve_tril) (TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) L,
                             const TEMPLATE(T, mat_t) B, int unit,
                             const TEMPLATE(T, ctx_t) ctx)
{
    if (B->r < TEMPLATE(CAP_T, MAT_SOLVE_TRI_ROWS_CUTOFF) ||
        B->c < TEMPLATE(CAP_T, MAT_SOLVE_TRI_COLS_CUTOFF))
    {
        TEMPLATE(T, mat_solve_tril_classical) (X, L, B, unit, ctx);
    }
    else
    {
        TEMPLATE(T, mat_solve_tril_recursive) (X, L, B, unit, ctx);
    }
}


#endif
