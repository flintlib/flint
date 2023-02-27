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
_TEMPLATE(T, poly_mul_classical) (TEMPLATE(T, struct) * rop,
                                  const TEMPLATE(T, struct) * op1, slong len1,
                                  const TEMPLATE(T, struct) * op2, slong len2,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    if (len1 == 1 && len2 == 1)
    {
        TEMPLATE(T, mul) (rop, op1, op2, ctx);
    }
    else
    {
        slong i;

        /* Set res[i] = poly1[i]*poly2[0] */
        _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (rop, op1, len1, op2, ctx);

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        _TEMPLATE(T, TEMPLATE(poly_scalar_mul, T)) (rop + len1, op2 + 1,
                                                    len2 - 1, op1 + len1 - 1,
                                                    ctx);

        /* out[i+j] += in1[i]*in2[j] */
        for (i = 0; i < len1 - 1; i++)
            _TEMPLATE(T, TEMPLATE(poly_scalar_addmul, T)) (rop + i + 1,
                                                           op2 + 1, len2 - 1,
                                                           op1 + i, ctx);
    }
}

void
TEMPLATE(T, poly_mul_classical) (TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op1,
                                 const TEMPLATE(T, poly_t) op2,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (rop == op1 || rop == op2)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, len, ctx);
        _TEMPLATE(T, poly_mul_classical) (t->coeffs, op1->coeffs, op1->length,
                                          op2->coeffs, op2->length, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, len, ctx);
        _TEMPLATE(T, poly_mul_classical) (rop->coeffs, op1->coeffs,
                                          op1->length, op2->coeffs,
                                          op2->length, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, len, ctx);
}


#endif
