/*
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
TEMPLATE(T, TEMPLATE(poly_set, T)) (TEMPLATE(T, poly_t) poly,
                                    const TEMPLATE(T, t) c,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    if (TEMPLATE(T, is_zero) (c, ctx))
    {
        TEMPLATE(T, poly_zero) (poly, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (poly, 1, ctx);
        TEMPLATE(T, set) (poly->coeffs + 0, c, ctx);
        _TEMPLATE(T, poly_set_length) (poly, 1, ctx);
    }
}


#endif
