/*
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
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
TEMPLATE(T, poly_init) (TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx)
{
    poly->coeffs = NULL;
    poly->alloc = 0;
    poly->length = 0;
}

void
TEMPLATE(T, poly_init2) (TEMPLATE(T, poly_t) poly, slong alloc,
                         const TEMPLATE(T, ctx_t) ctx)
{
    poly->coeffs = (alloc) ? _TEMPLATE(T, vec_init) (alloc, ctx) : NULL;
    poly->alloc = alloc;
    poly->length = 0;
}


#endif
