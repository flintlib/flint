/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart
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
_TEMPLATE(T, poly_div_basecase) (TEMPLATE(T, struct) * Q,
                                 TEMPLATE(T, struct) * R,
                                 const TEMPLATE(T, struct) * A, slong lenA,
                                 const TEMPLATE(T, struct) * B, slong lenB,
                                 const TEMPLATE(T, t) invB,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    const slong alloc = (R == NULL) ? lenA : 0;
    slong lenR = lenB - 1, iQ;

    if (alloc)
        R = _TEMPLATE(T, vec_init) (alloc, ctx);
    if (R != A)
        _TEMPLATE(T, vec_set) (R + lenR, A + lenR, lenA - lenR, ctx);

    for (iQ = lenA - lenB; iQ >= 0; iQ--)
    {
        if (TEMPLATE(T, is_zero) (R + lenA - 1, ctx))
        {
            TEMPLATE(T, zero) (Q + iQ, ctx);
        }
        else
        {
            TEMPLATE(T, mul) (Q + iQ, R + lenA - 1, invB, ctx);

            _TEMPLATE(T, TEMPLATE(vec_scalar_submul, T)) (R + lenA - lenR - 1,
                                                          B, lenR, Q + iQ,
                                                          ctx);
        }

        if (lenR - 1 >= iQ)
        {
            B++;
            lenR--;
        }

        lenA--;
    }

    if (alloc)
        _TEMPLATE(T, vec_clear) (R, alloc, ctx);
}

void
TEMPLATE(T, poly_div_basecase) (TEMPLATE(T, poly_t) Q,
                                const TEMPLATE(T, poly_t) A,
                                const TEMPLATE(T, poly_t) B,
                                const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    TEMPLATE(T, struct) * q;
    TEMPLATE(T, t) invB;

    if (lenA < lenB)
    {
        TEMPLATE(T, poly_zero) (Q, ctx);
        return;
    }

    TEMPLATE(T, init) (invB, ctx);
    TEMPLATE(T, inv) (invB, B->coeffs + (lenB - 1), ctx);

    if (Q == A || Q == B)
    {
        q = _TEMPLATE(T, vec_init) (lenQ, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (Q, lenQ, ctx);
        q = Q->coeffs;
    }

    _TEMPLATE(T, poly_div_basecase) (q, NULL, A->coeffs, lenA,
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

    TEMPLATE(T, clear) (invB, ctx);
}


#endif
