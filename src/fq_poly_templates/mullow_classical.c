/*
    Copyright (C) 2008, 2009 William Hart
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
_TEMPLATE(T, poly_mullow_classical) (
    TEMPLATE(T, struct) * rop,
    const TEMPLATE(T, struct) * op1, slong len1,
    const TEMPLATE(T, struct) * op2, slong len2,
    slong n, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_mullow_classical(rop, op1, len1, op2, len2, n, gr_ctx));
}

void
TEMPLATE(T, poly_mullow_classical) (TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op1,
                                    const TEMPLATE(T, poly_t) op2, slong n,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0 || n == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (n > len)
        n = len;

    if (rop == op1 || rop == op2)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, n, ctx);
        _TEMPLATE(T, poly_mullow_classical) (t->coeffs, op1->coeffs,
                                             op1->length, op2->coeffs,
                                             op2->length, n, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, n, ctx);
        _TEMPLATE(T, poly_mullow_classical) (rop->coeffs, op1->coeffs,
                                             op1->length, op2->coeffs,
                                             op2->length, n, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, n, ctx);
    _TEMPLATE(T, poly_normalise) (rop, ctx);
}


#endif
