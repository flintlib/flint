/*
    Copyright (C) 2008, 2009 William Hart
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
_TEMPLATE(T, TEMPLATE(poly_scalar_div, T)) (TEMPLATE(T, struct) * rop,
                                            const TEMPLATE(T, struct) * op,
                                            slong len, const TEMPLATE(T, t) x,
                                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        TEMPLATE(T, div) (rop + i, op + i, x, ctx);
}

void
TEMPLATE(T, TEMPLATE(poly_scalar_div, T)) (TEMPLATE(T, poly_t) rop,
                                           const TEMPLATE(T, poly_t) op,
                                           const TEMPLATE(T, t) x,
                                           const TEMPLATE(T, ctx_t) ctx)
{
    if (TEMPLATE(T, is_zero) (x, ctx))
    {
       flint_throw(FLINT_ERROR, "Exception (fq_poly_scalar_div) Division by zero");
    }
    if (TEMPLATE(T, poly_is_zero) (op, ctx))
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, op->length, ctx);
        _TEMPLATE(T, TEMPLATE(poly_scalar_div, T)) (rop->coeffs, op->coeffs,
                                                    op->length, x, ctx);
        _TEMPLATE(T, poly_set_length) (rop, op->length, ctx);
    }
}


#endif
