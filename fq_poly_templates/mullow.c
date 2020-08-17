/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
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
_TEMPLATE(T, poly_mullow) (TEMPLATE(T, struct) * rop,
                           const TEMPLATE(T, struct) * op1, slong len1,
                           const TEMPLATE(T, struct) * op2, slong len2,
                           slong n, const TEMPLATE(T, ctx_t) ctx)
{
    if (n < TEMPLATE(CAP_T, MULLOW_CLASSICAL_CUTOFF)
        || FLINT_MAX(len1, len2) < 6)
    {
        _TEMPLATE(T, poly_mullow_classical) (rop, op1, len1, op2, len2, n,
                                             ctx);
    }
    else
    {
        _TEMPLATE(T, poly_mullow_KS) (rop, op1, len1, op2, len2, n, ctx);
    }
}

void
TEMPLATE(T, poly_mullow) (TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op1,
                          const TEMPLATE(T, poly_t) op2,
                          slong n, const TEMPLATE(T, ctx_t) ctx)
{
    const slong len1 = op1->length;
    const slong len2 = op2->length;
    const slong lenr = op1->length + op2->length - 1;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
        return;
    }

    if (n > lenr)
        n = lenr;

    if (rop == op1 || rop == op2)
    {
        TEMPLATE(T, poly_t) t;

        TEMPLATE(T, poly_init2) (t, n, ctx);
        _TEMPLATE(T, poly_mullow) (t->coeffs, op1->coeffs, op1->length,
                                   op2->coeffs, op2->length, n, ctx);
        TEMPLATE(T, poly_swap) (rop, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, n, ctx);
        _TEMPLATE(T, poly_mullow) (rop->coeffs, op1->coeffs, op1->length,
                                   op2->coeffs, op2->length, n, ctx);
    }

    _TEMPLATE(T, poly_set_length) (rop, n, ctx);
    _TEMPLATE(T, poly_normalise) (rop, ctx);
}


#endif
