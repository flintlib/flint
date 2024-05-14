/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "long_extras.h"
#include "templates.h"

void
TEMPLATE(T, mat_init)(TEMPLATE(T, mat_t) mat, slong rows, slong cols,
                       const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    mat->r = rows;
    mat->c = cols;

    if (rows != 0)
        mat->rows = flint_malloc(rows * sizeof(TEMPLATE(T, struct) *));
    else
        mat->rows = NULL;

    if (rows != 0 && cols != 0)
    {
        slong j;
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = flint_malloc(num * sizeof(TEMPLATE(T, struct)));

        for (i = 0; i < rows; i++)
        {
            mat->rows[i] = mat->entries + i * cols;
            for (j = 0; j < cols; j++)
                TEMPLATE(T, init) (mat->rows[i] + j, ctx);
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
}
#endif
