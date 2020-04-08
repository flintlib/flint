/*
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
TEMPLATE(T, mat_mul_vec) (TEMPLATE(T, struct) *y,
                      const TEMPLATE(T, mat_t) A,
                      const TEMPLATE(T, struct) *x, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for (i = 0; i < A->r; ++i)
        _TEMPLATE(T, vec_dot) (&y[i], A->rows[i], x, A->c, ctx);
}


#endif
