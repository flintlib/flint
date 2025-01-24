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
    mat->entries = NULL;
    mat->r = rows;
    mat->c = cols;
    mat->stride = cols;

    if (rows != 0 && cols != 0)
    {
        slong num;
        int of;

        of = z_mul_checked(&num, rows, cols);

        if (of)
            flint_throw(FLINT_ERROR, "Overflow creating a %wd x %wd object\n", rows, cols);

        mat->entries = _TEMPLATE(T, vec_init)(num, ctx);
    }
}
#endif
