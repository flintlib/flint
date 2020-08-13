/*
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
_TEMPLATE(T, poly_neg) (TEMPLATE(T, struct) * rop,
                        const TEMPLATE(T, struct) * op, slong len,
                        const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        TEMPLATE(T, neg) (rop + i, op + i, ctx);
}

void
TEMPLATE(T, poly_neg) (TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                       const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_fit_length) (rop, op->length, ctx);
    _TEMPLATE(T, poly_neg) (rop->coeffs, op->coeffs, op->length, ctx);
    _TEMPLATE(T, poly_set_length) (rop, op->length, ctx);
}


#endif
