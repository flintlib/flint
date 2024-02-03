/*
    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2008, 2009 William Hart
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
_TEMPLATE(T, poly_div) (TEMPLATE(T, struct) * Q,
                                 const TEMPLATE(T, struct) * A, slong lenA,
                                 const TEMPLATE(T, struct) * B, slong lenB,
                                 const TEMPLATE(T, t) invB,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    if (lenB <= 15 || lenA - lenB <= 15)
        GR_MUST_SUCCEED(_gr_poly_div_basecase_preinv1(Q, A, lenA, B, lenB, invB, gr_ctx));
    else
        GR_MUST_SUCCEED(_gr_poly_div_newton(Q, A, lenA, B, lenB, gr_ctx));  /* todo: pass invB */
}

void
TEMPLATE(T, poly_div) (TEMPLATE(T, poly_t) Q,
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

    _TEMPLATE(T, poly_div) (q, A->coeffs, lenA,
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

/* flint 2.x compatibility needed by Nemo */
void
TEMPLATE(T, poly_div_basecase) (TEMPLATE(T, poly_t) Q,
                                const TEMPLATE(T, poly_t) A,
                                const TEMPLATE(T, poly_t) B,
                                const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_div) (Q, A, B, ctx);
}

#endif
