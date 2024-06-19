/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "templates.h"

void
TEMPLATE(T, mat_swap_entrywise)(TEMPLATE(T, mat_t) mat1,
		         TEMPLATE(T, mat_t) mat2, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;

    for (i = 0; i < TEMPLATE(T, mat_nrows)(mat1, ctx); i++)
        for (j = 0; j < TEMPLATE(T, mat_ncols)(mat1, ctx); j++)
            TEMPLATE(T, swap)(TEMPLATE(T, mat_entry)(mat2, i, j),
			      TEMPLATE(T, mat_entry)(mat1, i, j), ctx);
}

#endif
