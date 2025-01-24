/*
    Copyright (C) 2011 Fredrik Johansson
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
TEMPLATE(T, mat_neg) (TEMPLATE(T, mat_t) res,
                      const TEMPLATE(T, mat_t) mat,
                      const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (res->c < 1)
        return;

    for (i = 0; i < res->r; i++)
        _TEMPLATE(T, vec_neg) (TEMPLATE(T, mat_entry)(res, i, 0),
            TEMPLATE(T, mat_entry)(mat, i, 0), res->c, ctx);
}


#endif
