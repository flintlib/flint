/*
    Copyright (C) 2008, 2009, 2010 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
_TEMPLATE(T, vec_equal) (const TEMPLATE(T, struct) * vec1,
                         const TEMPLATE(T, struct) * vec2, slong len,
                         const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    if (vec1 == vec2)
        return 1;

    for (i = 0; i < len; i++)
        if (!TEMPLATE(T, equal) (vec1 + i, vec2 + i, ctx))
            return 0;

    return 1;
}


#endif
