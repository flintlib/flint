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
_TEMPLATE(T, poly_sqr) (TEMPLATE(T, struct) * rop,
                        const TEMPLATE(T, struct) * op, slong len,
                        const TEMPLATE(T, ctx_t) ctx)
{
    if (len < TEMPLATE(CAP_T, SQR_CLASSICAL_CUTOFF))
    {
        _TEMPLATE(T, poly_sqr_classical) (rop, op, len, ctx);
    }
#ifdef USE_SQR_REORDER
    else if (TEMPLATE(T, ctx_degree) (ctx) < 4)
    {
        _TEMPLATE(T, poly_sqr_reorder) (rop, op, len, ctx);
    }
#endif
    else
    {
        _TEMPLATE(T, poly_sqr_KS) (rop, op, len, ctx);
    }
}

void
TEMPLATE(T, poly_sqr) (TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                       const TEMPLATE(T, ctx_t) ctx)
{
    const slong rlen = 2 * op->length - 1;

    if (op->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (rop == op)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, rlen, ctx);
        _TEMPLATE(T, poly_sqr) (t->coeffs, op->coeffs, op->length, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, rlen, ctx);
        _TEMPLATE(T, poly_sqr) (rop->coeffs, op->coeffs, op->length, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, rlen, ctx);
}


#endif
