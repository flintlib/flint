/*
    Copyright (C) 2011 Sebastian Pancratz
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
_TEMPLATE(T, poly_divrem_f) (TEMPLATE(T, t) f,
                             TEMPLATE(T, struct) * Q, TEMPLATE(T, struct) * R,
                             const TEMPLATE(T, struct) * A, slong lenA,
                             const TEMPLATE(T, struct) * B, slong lenB,
                             const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) invB;

    TEMPLATE(T, init) (invB, ctx);
    TEMPLATE(T, gcdinv) (f, invB, B + lenB - 1, ctx);

    if (TEMPLATE(T, is_one) (f, ctx))
    {
        _TEMPLATE(T, poly_divrem) (Q, R, A, lenA, B, lenB, invB, ctx);
    }

    TEMPLATE(T, clear) (invB, ctx);
}

void
TEMPLATE(T, poly_divrem_f) (TEMPLATE(T, t) f,
                            TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                            const TEMPLATE(T, poly_t) A,
                            const TEMPLATE(T, poly_t) B,
                            const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;

    TEMPLATE(T, struct) * q, *r;
    TEMPLATE(T, t) invB;

    TEMPLATE(T, init) (invB, ctx);
    TEMPLATE(T, gcdinv) (f, invB, TEMPLATE(T, poly_lead) (B, ctx), ctx);

    if (!TEMPLATE(T, is_one) (f, ctx))
    {
        TEMPLATE(T, clear) (invB, ctx);
        return;
    }

    if (lenA < lenB)
    {
        TEMPLATE(T, poly_set) (R, A, ctx);
        TEMPLATE(T, poly_zero) (Q, ctx);
        TEMPLATE(T, clear) (invB, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        q = _TEMPLATE(T, vec_init) (lenQ, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (Q, lenQ, ctx);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _TEMPLATE(T, vec_init) (lenA, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (R, lenA, ctx);
        r = R->coeffs;
    }

    _TEMPLATE(T, poly_divrem_divconquer) (q, r, A->coeffs, lenA,
                                          B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _TEMPLATE(T, vec_clear) (Q->coeffs, Q->alloc, ctx);
        Q->coeffs = q;
        Q->alloc = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _TEMPLATE(T, poly_set_length) (Q, lenQ, ctx);
    }

    if (R == A || R == B)
    {
        _TEMPLATE(T, vec_clear) (R->coeffs, R->alloc, ctx);
        R->coeffs = r;
        R->alloc = lenA;
        R->length = lenA;
    }
    _TEMPLATE(T, poly_set_length) (R, lenB - 1, ctx);
    _TEMPLATE(T, poly_normalise) (R, ctx);

    TEMPLATE(T, clear) (invB, ctx);
}


#endif
