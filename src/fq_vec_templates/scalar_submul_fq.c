/*
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
_TEMPLATE3(T, vec_scalar_submul, T) (TEMPLATE(T, struct) * vec1,
                                     const TEMPLATE(T, struct) * vec2,
                                     slong len2, const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, t) y;

    TEMPLATE(T, init) (y, ctx);

    for (i = 0; i < len2; i++)
    {
        TEMPLATE(T, mul) (y, vec2 + i, x, ctx);
        TEMPLATE(T, sub) (vec1 + i, vec1 + i, y, ctx);
    }

    TEMPLATE(T, clear) (y, ctx);
}


#endif
