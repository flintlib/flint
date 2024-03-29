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

#include "templates.h"

void
_TEMPLATE(T, vec_clear) (TEMPLATE(T, struct) * vec, slong len,
                         const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        TEMPLATE(T, clear) (vec + i, ctx);
    flint_free(vec);
}

#endif
