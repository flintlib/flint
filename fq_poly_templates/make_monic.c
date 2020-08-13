/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart
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
_TEMPLATE(T, poly_make_monic) (TEMPLATE(T, struct) * rop,
                               const TEMPLATE(T, struct) * op, slong length,
                               const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) inv;
    TEMPLATE(T, init) (inv, ctx);
    TEMPLATE(T, inv) (inv, &op[length - 1], ctx);
    _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (rop, op, length, inv, ctx);
    TEMPLATE(T, clear) (inv, ctx);
}

void
TEMPLATE(T, poly_make_monic) (TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op,
                              const TEMPLATE(T, ctx_t) ctx)
{
    if (op->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    TEMPLATE(T, poly_fit_length) (rop, op->length, ctx);
    _TEMPLATE(T, poly_make_monic) (rop->coeffs, op->coeffs, op->length, ctx);
    _TEMPLATE(T, poly_set_length) (rop, op->length, ctx);
}


#endif
