/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, vec_scalar_mul_si) (TEMPLATE(T, struct) * vec1,
                                  const TEMPLATE(T, struct) * vec2,
                                  slong len2, const slong x,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < len2; i++)
    {
        TEMPLATE(T, mul_si) (vec1 + i, vec2 + i, x, ctx);
    }
}


#endif
