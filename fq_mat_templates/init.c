/*
    Copyright (C) 2010 William Hart
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
TEMPLATE(T, mat_init) (TEMPLATE(T, mat_t) mat, slong rows, slong cols,
                       const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (rows != 0)
        mat->rows = (TEMPLATE(T, struct) **) flint_malloc(rows
			                       * sizeof(TEMPLATE(T, struct) *));
    else
        mat->rows = NULL;

    if (rows != 0 && cols != 0)       /* Allocate space for r*c small entries */
    {
        slong j;

	mat->entries = (TEMPLATE(T, struct) *) flint_malloc(flint_mul_sizes(rows, cols)
                                                * sizeof(TEMPLATE(T, struct)));

        for (i = 0; i < rows; i++)
        {
            mat->rows[i] = mat->entries + i * cols;
            for (j = 0; j < cols; j++)
            {
                TEMPLATE(T, init) (mat->rows[i] + j, ctx);
            }
        }
    }
    else
    {
        mat->entries = NULL;
	if (rows != 0)
	{
            for (i = 0; i < rows; i++)
                mat->rows[i] = NULL;
	}
    }

    mat->r = rows;
    mat->c = cols;
}


#endif
