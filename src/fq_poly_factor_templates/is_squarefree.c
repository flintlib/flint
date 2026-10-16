/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
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

int
_TEMPLATE(T, poly_is_squarefree) (const TEMPLATE(T, struct) * f, slong len,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    gr_poly_struct t;
    truth_t res;

    if (len <= 2)
        return len != 0;

    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    t.coeffs = (TEMPLATE(T, struct) *) f;
    t.length = len;
    t.alloc = len;

    res = gr_poly_is_squarefree(&t, gr_ctx);

    if (res == T_UNKNOWN)
        flint_throw(FLINT_ERROR, "Exception in poly_is_squarefree: "
                                 "unable to decide\n");

    return res == T_TRUE;
}

int
TEMPLATE(T, poly_is_squarefree) (const TEMPLATE(T, poly_t) f,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_is_squarefree) (f->coeffs, f->length, ctx);
}


#endif
