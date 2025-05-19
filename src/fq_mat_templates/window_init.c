/*
    Copyright (C) 2008, 2009 William Hart.
    Copyright (C) 2008, Richard Howell-Peak
    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2010, Fredrik Johansson
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
TEMPLATE(T, mat_window_init) (TEMPLATE(T, mat_t) window,
                              const TEMPLATE(T, mat_t) mat,
                              slong r1, slong c1, slong r2, slong c2,
                              const TEMPLATE(T, ctx_t) FLINT_UNUSED(ctx))
{
    FLINT_ASSERT(r1 >= 0 && r1 <= r2 && r2 <= mat->r);
    FLINT_ASSERT(c2 >= 0 && c1 <= c2 && c2 <= mat->c);

    window->entries = TEMPLATE(T, mat_entry)(mat, r1, c1);
    window->r = r2 - r1;
    window->c = c2 - c1;
    window->stride = mat->stride;
}


#endif
