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

TEMPLATE(T, struct)*
_TEMPLATE(T, vec_init) (slong len, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    TEMPLATE(T, struct) * v;
    v = flint_malloc(len * sizeof(TEMPLATE(T, struct)));
    for (i = 0; i < len; i++)
        TEMPLATE(T, init) (v + i, ctx);
    return v;
}


#endif
