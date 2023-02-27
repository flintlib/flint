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
TEMPLATE(T, poly_set_coeff) (TEMPLATE(T, poly_t) poly, slong n,
                             const TEMPLATE(T, t) x,
                             const TEMPLATE(T, ctx_t) ctx)
{
    if (TEMPLATE(T, is_zero) (x, ctx))
    {
        if (n >= poly->length)
            return;

        TEMPLATE(T, zero) (poly->coeffs + n, ctx);

        if (n == poly->length - 1) /* only necessary when setting leading coefficient */
            _TEMPLATE(T, poly_normalise) (poly, ctx);
    }
    else
    {
        slong i;

        TEMPLATE(T, poly_fit_length) (poly, n + 1, ctx);

        if (n + 1 > poly->length)   /* Insert zeros if needed */
        {
            for (i = poly->length; i < n; i++)
                TEMPLATE(T, zero) (poly->coeffs + i, ctx);
            poly->length = n + 1;
        }
    }

    TEMPLATE(T, set) (poly->coeffs + n, x, ctx);
}


#endif
