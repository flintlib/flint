/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

ulong
TEMPLATE(T, poly_remove) (TEMPLATE(T, poly_t) f, const TEMPLATE(T, poly_t) g,
                          const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) q, r;
    ulong i = 0;

    TEMPLATE(T, poly_init) (q, ctx);
    TEMPLATE(T, poly_init) (r, ctx);

    while (1)
    {
        if (f->length < g->length)
            break;
        TEMPLATE(T, poly_divrem) (q, r, f, g, ctx);
        if (r->length == 0)
            TEMPLATE(T, poly_swap) (q, f, ctx);
        else
            break;
        i++;
    }

    TEMPLATE(T, poly_clear) (q, ctx);
    TEMPLATE(T, poly_clear) (r, ctx);

    return i;
}


#endif
