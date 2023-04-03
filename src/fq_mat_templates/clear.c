/*
    Copyright (C) 2010 William Hart
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
TEMPLATE(T, mat_clear) (TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx)
{
    if (mat->entries)
    {
        slong i;
        for (i = 0; i < mat->r * mat->c; i++)
            TEMPLATE(T, clear) (mat->entries + i, ctx); /* Clear all coefficients */
        flint_free(mat->entries);   /* Clean up array of entries */
        flint_free(mat->rows);  /* Clean up row array */
    } else if (mat->r != 0)
        flint_free(mat->rows);
}


#endif
