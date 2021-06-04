/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_scalar_mul) (TEMPLATE(T, mat_t) B,
                         const TEMPLATE(T, mat_t) A,
                         const TEMPLATE(T, t) c,
                         const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for (i = 0; i < B->r; i++)
    {
        _TEMPLATE(T, TEMPLATE(vec_scalar_mul, T)) (B->rows[i], A->rows[i], A->c, c, ctx);
    }
}


#endif
