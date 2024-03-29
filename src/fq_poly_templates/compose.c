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
_TEMPLATE(T, poly_compose) (TEMPLATE(T, struct) * rop,
                            const TEMPLATE(T, struct) * op1, slong len1,
                            const TEMPLATE(T, struct) * op2, slong len2,
                            const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_compose(rop, op1, len1, op2, len2, gr_ctx));
}

void
TEMPLATE(T, poly_compose) (TEMPLATE(T, poly_t) rop,
                           const TEMPLATE(T, poly_t) op1,
                           const TEMPLATE(T, poly_t) op2,
                           const TEMPLATE(T, ctx_t) ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong lenr = (len1 - 1) * (len2 - 1) + 1;

    if (len1 == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        TEMPLATE(T, TEMPLATE(poly_set, T)) (rop, op1->coeffs + 0, ctx);
    }
    else if (rop != op1 && rop != op2)
    {
        TEMPLATE(T, poly_fit_length) (rop, lenr, ctx);
        _TEMPLATE(T, poly_compose) (rop->coeffs, op1->coeffs, len1,
                                    op2->coeffs, len2, ctx);
        _TEMPLATE(T, poly_set_length) (rop, lenr, ctx);
        _TEMPLATE(T, poly_normalise) (rop, ctx);
    }
    else
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, lenr, ctx);
        _TEMPLATE(T, poly_compose) (t->coeffs, op1->coeffs, len1, op2->coeffs,
                                    len2, ctx);
        _TEMPLATE(T, poly_set_length) (t, lenr, ctx);
        _TEMPLATE(T, poly_normalise) (t, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
}


#endif
