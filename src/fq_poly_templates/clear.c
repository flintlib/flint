/*
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
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
TEMPLATE(T, poly_clear) (TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    if (poly->coeffs)
    {
        _TEMPLATE(T, vec_clear) (poly->coeffs, poly->alloc, ctx);
    }
}


#endif
