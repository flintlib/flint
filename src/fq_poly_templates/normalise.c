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
_TEMPLATE(T, poly_normalise) (TEMPLATE(T, poly_t) poly,
                              const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = poly->length - 1;
         (i >= 0) && TEMPLATE(T, is_zero) (poly->coeffs + i, ctx); i--) ;
    poly->length = i + 1;
}

void
_TEMPLATE(T, poly_normalise2) (const TEMPLATE(T, struct) * poly,
		slong * length, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    for (i = (*length) - 1; (i >= 0) && TEMPLATE(T, is_zero) (poly + i, ctx);
         i--) ;
    (*length) = i + 1;
}


#endif
