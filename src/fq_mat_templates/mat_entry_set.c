/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_entry_set)(TEMPLATE(T, mat_t) mat, slong i, slong j,
                           const TEMPLATE(T, t) x,
                           const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, set)(TEMPLATE(T, mat_entry)(mat, i, j), x, ctx);
}

#endif
