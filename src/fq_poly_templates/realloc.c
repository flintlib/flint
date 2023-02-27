/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_realloc) (TEMPLATE(T, poly_t) poly, slong alloc,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        TEMPLATE(T, poly_clear) (poly, ctx);
        TEMPLATE(T, poly_init) (poly, ctx);
    }
    else if (poly->alloc)       /* Realloc */
    {
        for (i = alloc; i < poly->alloc; i++)
            TEMPLATE(T, clear) (poly->coeffs + i, ctx);

        poly->coeffs =
            (TEMPLATE(T, struct) *) flint_realloc(poly->coeffs,
                                                  alloc *
                                                  sizeof(TEMPLATE(T, struct)));

        for (i = poly->alloc; i < alloc; i++)
            TEMPLATE(T, init) (poly->coeffs + i, ctx);

        poly->length = FLINT_MIN(poly->length, alloc);
        _TEMPLATE(T, poly_normalise) (poly, ctx);
    }
    else                        /* Nothing allocated already so do it now */
    {
        poly->coeffs =
            (TEMPLATE(T, struct) *) flint_malloc(alloc *
                                                 sizeof(TEMPLATE(T, struct)));

        for (i = 0; i < alloc; i++)
            TEMPLATE(T, init) (poly->coeffs + i, ctx);
    }
    poly->alloc = alloc;
}


#endif
