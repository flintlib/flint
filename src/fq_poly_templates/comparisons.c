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

int
TEMPLATE(T, poly_is_one)(const TEMPLATE(T, poly_t) op,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return (op->length == 1) && (TEMPLATE(T, is_one)(op->coeffs + 0, ctx));
}

int
TEMPLATE(T, poly_is_unit)(const TEMPLATE(T, poly_t) op,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return (op->length == 1) && (!(TEMPLATE(T, is_zero)(op->coeffs + 0, ctx)));
}

int
TEMPLATE3(T, poly_equal, T)(const TEMPLATE(T, poly_t) poly,
                            const TEMPLATE(T, t) c,
                            const TEMPLATE(T, ctx_t) ctx)
{
    return ((poly->length == 0) && TEMPLATE(T, is_zero)(c, ctx)) ||
        ((poly->length == 1) && TEMPLATE(T, equal)(poly->coeffs, c, ctx));
}

#endif
