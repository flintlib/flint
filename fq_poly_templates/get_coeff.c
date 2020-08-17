/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
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
TEMPLATE(T, poly_get_coeff) (TEMPLATE(T, t) x, const TEMPLATE(T, poly_t) poly,
                             slong n, const TEMPLATE(T, ctx_t) ctx)
{
    if (n < poly->length)
        TEMPLATE(T, set) (x, poly->coeffs + n, ctx);
    else
        TEMPLATE(T, zero) (x, ctx);
}


#endif
