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

/* Assumes op1 and op2 are not length 0. */
void
_TEMPLATE(T, poly_mulhigh_classical) (
    TEMPLATE(T, struct) * rop,
    const TEMPLATE(T, struct) * op1, slong len1,
    const TEMPLATE(T, struct) * op2, slong len2,
    slong start, const TEMPLATE(T, ctx_t) ctx)
{
    slong m, n;

    _TEMPLATE(T, vec_zero) (rop, start, ctx);

    if (len1 == 1)              /* Special case if the length of both inputs is 1 */
    {
        if (start == 0)
            TEMPLATE(T, mul) (rop, op1, op2, ctx);
    }
    else                        /* Ordinary case */
    {
        slong i;

        /* Set res[i] = poly1[i]*poly2[0] */
        if (start < len1)
            _TEMPLATE3(T, vec_scalar_mul, T) (rop + start, op1 + start,
                                              len1 - start, op2, ctx);

        if (len2 == 1)
            return;

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        m = FLINT_MAX(len1 - 1, start);
        _TEMPLATE3(T, vec_scalar_mul, T) (rop + m, op2 + m - len1 + 1,
                                          len2 - 1 + len1 - m, op1 + len1 - 1,
                                          ctx);

        /* out[i+j] += in1[i]*in2[j] */
        m = FLINT_MAX(start, len2 - 1);
        for (i = m - len2 + 1; i < len1 - 1; i++)
        {
            n = FLINT_MAX(i + 1, start);
            _TEMPLATE3(T, vec_scalar_addmul, T) (rop + n, op2 + n - i,
                                                 len2 + i - n, op1 + i, ctx);
        }
    }
}

void
TEMPLATE(T, poly_mulhigh_classical) (TEMPLATE(T, poly_t) rop,
                                     const TEMPLATE(T, poly_t) op1,
                                     const TEMPLATE(T, poly_t) op2,
                                     slong start, const TEMPLATE(T, ctx_t) ctx)
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
