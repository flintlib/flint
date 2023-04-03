/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_set_fmpz_mod_mat) (TEMPLATE(T, mat_t) mat1,
		       const fmpz_mod_mat_t mat2, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    TEMPLATE(T, t) t;

    TEMPLATE(T, init)(t, ctx);

    for (i = 0; i < mat1->r; i++)
    {
        for (j = 0; j < mat1->c; j++)
        {
            TEMPLATE(T, set_fmpz)(t, fmpz_mod_mat_entry(mat2, i, j), ctx);
            TEMPLATE(T, mat_entry_set)(mat1, i, j, t, ctx);
        }
    }

    TEMPLATE(T, clear)(t, ctx);
}

#endif
