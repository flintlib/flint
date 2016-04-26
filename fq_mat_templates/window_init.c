/*
    Copyright (C) 2008, 2009 William Hart.
    Copyright (C) 2008, Richard Howell-Peak
    Copyright (C) 2008, Martin Albrecht
    Copyright (C) 2010, Fredrik Johansson
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
TEMPLATE(T, mat_window_init) (TEMPLATE(T, mat_t) window,
                              const TEMPLATE(T, mat_t) mat,
                              slong r1, slong c1, slong r2, slong c2,
                              const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    window->entries = NULL;

    window->rows = flint_malloc((r2 - r1) * sizeof(TEMPLATE(T, struct) *));

    if (mat->c > 0)
    {
        for (i = 0; i < r2 - r1; i++)
            window->rows[i] = mat->rows[r1 + i] + c1;
    }

    window->r = r2 - r1;
    window->c = c2 - c1;
}


#endif
