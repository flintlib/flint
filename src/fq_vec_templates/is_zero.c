/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
_TEMPLATE(T, vec_is_zero) (const TEMPLATE(T, struct) * vec, slong len,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!TEMPLATE(T, is_zero) (vec + i, ctx))
            return 0;
    return 1;
}


#endif
