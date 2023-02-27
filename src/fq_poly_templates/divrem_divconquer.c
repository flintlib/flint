/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

static void
__TEMPLATE(T, poly_divrem_divconquer) (TEMPLATE(T, struct) * Q,
                                       TEMPLATE(T, struct) * R,
                                       const TEMPLATE(T, struct) * A, slong lenA,
                                       const TEMPLATE(T, struct) * B, slong lenB,
                                       const TEMPLATE(T, t) invB,
                                       const TEMPLATE(T, ctx_t) ctx)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const slong n1 = lenA - lenB + 1;
        const slong n2 = lenB - n1;

        const TEMPLATE(T, struct) * p1 = A + n2;
        const TEMPLATE(T, struct) * d1 = B + n2;
        const TEMPLATE(T, struct) * d2 = B;

        TEMPLATE(T, struct) * W =
            _TEMPLATE(T, vec_init) ((2 * n1 - 1) + lenB - 1, ctx);

        TEMPLATE(T, struct) * d1q1 = R + n2;
        TEMPLATE(T, struct) * d2q1 = W + (2 * n1 - 1);

        _TEMPLATE(T, poly_divrem_divconquer_recursive) (Q, d1q1, W, p1, d1, n1,
                                                        invB, ctx);

        /*
           Compute d2q1 = Q d2, of length lenB - 1
         */

        if (n1 >= n2)
            _TEMPLATE(T, poly_mul) (d2q1, Q, n1, d2, n2, ctx);
        else
            _TEMPLATE(T, poly_mul) (d2q1, d2, n2, Q, n1, ctx);

        /*
           Compute BQ = d1q1 * x^n1 + d2q1, of length lenB - 1; 
           then compute R = A - BQ
         */

        _TEMPLATE(T, vec_swap) (R, d2q1, n2, ctx);
        _TEMPLATE(T, poly_add) (R + n2, R + n2, n1 - 1, d2q1 + n2, n1 - 1,
                                ctx);
        _TEMPLATE(T, poly_sub) (R, A, lenA, R, lenA, ctx);

        _TEMPLATE(T, vec_clear) (W, (2 * n1 - 1) + lenB - 1, ctx);
    }
    else                        /* lenA = 2 * lenB - 1 */
    {
        TEMPLATE(T, struct) * W = _TEMPLATE(T, vec_init) (lenA, ctx);

        _TEMPLATE(T, poly_divrem_divconquer_recursive) (Q, R, W, A, B, lenB,
                                                        invB, ctx);

        _TEMPLATE(T, poly_sub) (R, A, lenB - 1, R, lenB - 1, ctx);

        _TEMPLATE(T, vec_clear) (W, lenA, ctx);
    }
}

void
_TEMPLATE(T, poly_divrem_divconquer) (TEMPLATE(T, struct) * Q,
                                      TEMPLATE(T, struct) * R,
                                      const TEMPLATE(T, struct) * A,
                                      slong lenA, const TEMPLATE(T,
                                                                 struct) * B,
                                      slong lenB, const TEMPLATE(T, t) invB,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    if (lenA <= 2 * lenB - 1)
    {
        __TEMPLATE(T, poly_divrem_divconquer) (Q, R, A, lenA, B, lenB, invB,
                                               ctx);
    }
    else                        /* lenA > 2 * lenB - 1 */
    {
        slong shift, n = 2 * lenB - 1;
        TEMPLATE(T, struct) * QB, *W;

        _TEMPLATE(T, vec_set) (R, A, lenA, ctx);
        W = _TEMPLATE(T, vec_init) (2 * n, ctx);
        QB = W + n;

        while (lenA >= n)
        {
            shift = lenA - n;
            _TEMPLATE(T, poly_divrem_divconquer_recursive) (Q + shift, QB,
                                                            W, R + shift, B,
                                                            lenB, invB, ctx);
            _TEMPLATE(T, poly_sub) (R + shift, R + shift, n, QB, n, ctx);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            __TEMPLATE(T, poly_divrem_divconquer) (Q, W, R, lenA, B, lenB,
                                                   invB, ctx);
            _TEMPLATE(T, vec_swap) (W, R, lenA, ctx);
        }

        _TEMPLATE(T, vec_clear) (W, 2 * n, ctx);
    }
}

void
TEMPLATE(T, poly_divrem_divconquer) (TEMPLATE(T, poly_t) Q,
                                     TEMPLATE(T, poly_t) R,
                                     const TEMPLATE(T, poly_t) A,
                                     const TEMPLATE(T, poly_t) B,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    const slong lenA = A->length;
    const slong lenB = B->length;
    const slong lenQ = lenA - lenB + 1;

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
