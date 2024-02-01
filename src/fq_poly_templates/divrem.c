/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Fredrik Johansson

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
_TEMPLATE(T, poly_divrem) (TEMPLATE(T, struct) * Q,
                                    TEMPLATE(T, struct) * R,
                                    const TEMPLATE(T, struct) * A, slong lenA,
                                    const TEMPLATE(T, struct) * B, slong lenB,
                                    const TEMPLATE(T, t) invB,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    /* todo: tuning, and make sure it is also used in gr_poly_divrem */
    if (lenB <= 10 || lenA - lenB <= 1)
        GR_MUST_SUCCEED(_gr_poly_divrem_basecase_preinv1(Q, R, A, lenA, B, lenB, invB, gr_ctx));
    else
        GR_MUST_SUCCEED(_gr_poly_divrem_newton(Q, R, A, lenA, B, lenB, gr_ctx));  /* todo: pass invB */
}

void
TEMPLATE(T, poly_divrem) (TEMPLATE(T, poly_t) Q,
                                   TEMPLATE(T, poly_t) R,
                                   const TEMPLATE(T, poly_t) A,
                                   const TEMPLATE(T, poly_t) B,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    TEMPLATE(T, struct) * q, *r;
    TEMPLATE(T, t) invB;

    if (lenA < lenB)
    {
        TEMPLATE(T, poly_set) (R, A, ctx);
        TEMPLATE(T, poly_zero) (Q, ctx);
        return;
    }

    TEMPLATE(T, init) (invB, ctx);
    TEMPLATE(T, inv) (invB, TEMPLATE(T, poly_lead) (B, ctx), ctx);

    if (Q == A || Q == B)
    {
        q = _TEMPLATE(T, vec_init) (lenQ, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (Q, lenQ, ctx);
        q = Q->coeffs;
    }
    if (R == B)
    {
        r = _TEMPLATE(T, vec_init) (lenA, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (R, lenA, ctx);
        r = R->coeffs;
    }

    _TEMPLATE(T, poly_divrem) (q, r, A->coeffs, lenA,
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
    if (R == B)
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
