/*
    Copyright (C) 2010 William Hart
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
_TEMPLATE(T, poly_mulhigh) (TEMPLATE(T, struct) * rop,
                            const TEMPLATE(T, struct) * op1, slong len1,
                            const TEMPLATE(T, struct) * op2, slong len2,
                            slong n, TEMPLATE(T, ctx_t) ctx)
{
    if (FLINT_MAX(len1, len2) < 6)
    {
        _TEMPLATE(T, poly_mulhigh_classical) (rop, op1, len1, op2, len2, n,
                                              ctx);
    }
    else
    {
        _TEMPLATE(T, poly_mul_KS) (rop, op1, len1, op2, len2, ctx);
    }
}

void
TEMPLATE(T, poly_mulhigh) (TEMPLATE(T, poly_t) rop,
                           const TEMPLATE(T, poly_t) op1,
                           const TEMPLATE(T, poly_t) op2, slong start,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong len_out = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0 || start >= len_out)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        TEMPLATE(T, poly_t) temp;
        TEMPLATE(T, poly_init2) (temp, len_out, ctx);
        if (op1->length >= op2->length)
            _TEMPLATE(T, poly_mulhigh_classical) (temp->coeffs, op1->coeffs,
                                                  op1->length, op2->coeffs,
                                                  op2->length, start, ctx);
        else
            _TEMPLATE(T, poly_mulhigh_classical) (temp->coeffs, op2->coeffs,
                                                  op2->length, op1->coeffs,
                                                  op1->length, start, ctx);
        TEMPLATE(T, poly_swap) (rop, temp, ctx);
        TEMPLATE(T, poly_clear) (temp, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, len_out, ctx);
        if (op1->length >= op2->length)
            _TEMPLATE(T, poly_mulhigh_classical) (rop->coeffs, op1->coeffs,
                                                  op1->length, op2->coeffs,
                                                  op2->length, start, ctx);
        else
            _TEMPLATE(T, poly_mulhigh_classical) (rop->coeffs, op2->coeffs,
                                                  op2->length, op1->coeffs,
                                                  op1->length, start, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, len_out, ctx);
    _TEMPLATE(T, poly_normalise) (rop, ctx);
}


#endif
