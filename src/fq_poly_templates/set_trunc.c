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
TEMPLATE(T, poly_set_trunc) (TEMPLATE(T, poly_t) poly1, TEMPLATE(T, poly_t) poly2, slong len,
                            const TEMPLATE(T, ctx_t) ctx)
{
    if (poly1 == poly2)
        TEMPLATE(T, poly_truncate) (poly1, len, ctx);
    else if (len >= poly2->length)
        TEMPLATE(T, poly_set) (poly1, poly2, ctx);
    else
    {
        slong i;

        TEMPLATE(T, poly_fit_length) (poly1, len, ctx);

        for (i = 0; i < len; i++)
            TEMPLATE(T, set) (poly1->coeffs + i, poly2->coeffs + i, ctx);
        _TEMPLATE(T, poly_set_length) (poly1, len, ctx);
        _TEMPLATE(T, poly_normalise) (poly1, ctx);
    }
}


#endif
