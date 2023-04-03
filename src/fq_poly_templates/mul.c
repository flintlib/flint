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
_TEMPLATE(T, poly_mul) (TEMPLATE(T, struct) * rop,
                        const TEMPLATE(T, struct) * op1, slong len1,
                        const TEMPLATE(T, struct) * op2, slong len2,
                        const TEMPLATE(T, ctx_t) ctx)
{
    if (FLINT_MAX(len1, len2) < TEMPLATE(CAP_T, MUL_CLASSICAL_CUTOFF))
    {
        _TEMPLATE(T, poly_mul_classical) (rop, op1, len1, op2, len2, ctx);
    }
#ifdef USE_MUL_REORDER
    else if (TEMPLATE(T, ctx_degree) (ctx) < 4)
    {
        _TEMPLATE(T, poly_mul_reorder) (rop, op1, len1, op2, len2, ctx);
    }
#endif
    else
    {
        _TEMPLATE(T, poly_mul_KS) (rop, op1, len1, op2, len2, ctx);
    }
}

void
TEMPLATE(T, poly_mul) (TEMPLATE(T, poly_t) rop,
                       const TEMPLATE(T, poly_t) op1,
                       const TEMPLATE(T, poly_t) op2,
                       const TEMPLATE(T, ctx_t) ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong rlen = op1->length + op2->length - 1;

    if (len1 == 0 || len2 == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, rlen, ctx);
        _TEMPLATE(T, poly_mul) (t->coeffs, op1->coeffs, len1, op2->coeffs,
                                len2, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, rlen, ctx);
        _TEMPLATE(T, poly_mul) (rop->coeffs, op1->coeffs, len1, op2->coeffs,
                                len2, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, rlen, ctx);
}


#endif
