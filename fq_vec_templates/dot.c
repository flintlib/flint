/*
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
_TEMPLATE(T, vec_dot) (TEMPLATE(T, t) res,
                       const TEMPLATE(T, struct) * vec1,
                       const TEMPLATE(T, struct) * vec2,
                       slong len2, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, t) x;
    TEMPLATE(T, init) (x, ctx);

    TEMPLATE(T, zero) (res, ctx);

    for (i = 0; i < len2; i++)
    {
        TEMPLATE(T, mul) (x, vec1 + i, vec2 + i, ctx);
        TEMPLATE(T, add) (res, res, x, ctx);
    }

    TEMPLATE(T, clear) (x, ctx);
}

#endif
