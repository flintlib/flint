/*
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
#include "gr_poly.h"

void
_TEMPLATE3(T, poly_evaluate, T) (TEMPLATE(T, t) rop,
                                 const TEMPLATE(T, struct) * op, slong len,
                                 const TEMPLATE(T, t) a,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_evaluate(rop, op, len, a, gr_ctx));
}

void
TEMPLATE3(T, poly_evaluate, T) (TEMPLATE(T, t) rop,
                                const TEMPLATE(T, poly_t) f,
                                const TEMPLATE(T, t) a,
                                const TEMPLATE(T, ctx_t) ctx)
{
    if (rop == a)
    {
        TEMPLATE(T, t) t;
        TEMPLATE(T, init) (t, ctx);
        _TEMPLATE3(T, poly_evaluate, T) (t, f->coeffs, f->length, a, ctx);
        TEMPLATE(T, swap) (rop, t, ctx);
        TEMPLATE(T, clear) (t, ctx);
    }
    else
    {
        _TEMPLATE3(T, poly_evaluate, T) (rop, f->coeffs, f->length, a, ctx);
    }
}


#endif
