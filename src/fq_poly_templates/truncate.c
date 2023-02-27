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
TEMPLATE(T, poly_truncate) (TEMPLATE(T, poly_t) poly, slong len,
                            const TEMPLATE(T, ctx_t) ctx)
{
    if (poly->length > len)
    {
        slong i;

        for (i = len; i < poly->length; i++)
            TEMPLATE(T, zero) (poly->coeffs + i, ctx);
        poly->length = len;
        _TEMPLATE(T, poly_normalise) (poly, ctx);
    }
}


#endif
