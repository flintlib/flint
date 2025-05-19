/*
    Copyright (C) 2008-2009 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_set) (TEMPLATE(T, mat_t) mat1, const TEMPLATE(T, mat_t) mat2,
                      const TEMPLATE(T, ctx_t) ctx)
{
    if (mat1 != mat2)
    {
        slong i;

        if (mat2->r && mat2->c)
            for (i = 0; i < mat2->r; i++)
                _TEMPLATE(T, vec_set) (TEMPLATE(T, mat_entry)(mat1, i, 0),
                    TEMPLATE(T, mat_entry)(mat2, i, 0), mat2->c, ctx);
    }
}


#endif
