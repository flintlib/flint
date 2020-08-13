/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdlib.h>

slong
_TEMPLATE(T, poly_gcd_euclidean) (TEMPLATE(T, struct) * G,
                                  const TEMPLATE(T, struct) * A, slong lenA,
                                  const TEMPLATE(T, struct) * B, slong lenB,
                                  const TEMPLATE(T, t) invB,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    if (lenB == 1)
    {
        TEMPLATE(T, one) (G, ctx);
        return 1;
    }
    else                        /* lenA >= lenB > 1 */
    {
        const slong lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
        TEMPLATE(T, t) invR3;
        TEMPLATE(T, struct) * Q, *R1, *R2, *R3, *T, *W;
        slong lenR2, lenR3;

        W = _TEMPLATE(T, vec_init) (lenW, ctx);
        Q = W;
        R1 = W + FLINT_MAX(lenA - lenB + 1, lenB);
        R2 = R1 + lenA;
        R3 = R2 + lenB;

        _TEMPLATE(T, poly_divrem) (Q, R1, A, lenA, B, lenB, invB, ctx);

        lenR3 = lenB - 1;
        TEMPLATE(CAP_T, VEC_NORM) (R1, lenR3, ctx);

        if (lenR3 == 0)
        {
            _TEMPLATE(T, vec_set) (G, B, lenB, ctx);
            _TEMPLATE(T, vec_clear) (W, lenW, ctx);
            return lenB;
        }

        TEMPLATE(T, init) (invR3, ctx);

        T = R3;
        R3 = R1;
        R1 = T;
        _TEMPLATE(T, vec_set) (R2, B, lenB, ctx);
        lenR2 = lenB;

        do
        {
            TEMPLATE(T, inv) (invR3, R3 + (lenR3 - 1), ctx);

            _TEMPLATE(T, poly_divrem) (Q, R1, R2, lenR2, R3, lenR3, invR3,
                                       ctx);
            lenR2 = lenR3--;
            TEMPLATE(CAP_T, VEC_NORM) (R1, lenR3, ctx);
            T = R2;
            R2 = R3;
            R3 = R1;
            R1 = T;
        }
        while (lenR3 > 0);

        _TEMPLATE(T, vec_set) (G, R2, lenR2, ctx);

        _TEMPLATE(T, vec_clear) (W, lenW, ctx);
        TEMPLATE(T, clear) (invR3, ctx);

        return lenR2;
    }

}

void
TEMPLATE(T, poly_gcd_euclidean) (TEMPLATE(T, poly_t) G,
                                 const TEMPLATE(T, poly_t) A,
                                 const TEMPLATE(T, poly_t) B,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    if (A->length < B->length)
    {
        TEMPLATE(T, poly_gcd_euclidean) (G, B, A, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        TEMPLATE(T, t) invB;
        TEMPLATE(T, struct) * g;

        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            TEMPLATE(T, poly_zero) (G, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            TEMPLATE(T, poly_make_monic) (G, A, ctx);
        }
        else                    /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                g = _TEMPLATE(T, vec_init) (FLINT_MIN(lenA, lenB), ctx);
            }
            else
            {
                TEMPLATE(T, poly_fit_length) (G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            TEMPLATE(T, init) (invB, ctx);
            TEMPLATE(T, inv) (invB, TEMPLATE(T, poly_lead) (B, ctx), ctx);
            lenG = _TEMPLATE(T, poly_gcd_euclidean) (g, A->coeffs, lenA,
                                                     B->coeffs, lenB, invB,
                                                     ctx);
            TEMPLATE(T, clear) (invB, ctx);

            if (G == A || G == B)
            {
                _TEMPLATE(T, vec_clear) (G->coeffs, G->alloc, ctx);
                G->coeffs = g;
                G->alloc = FLINT_MIN(lenA, lenB);
                G->length = FLINT_MIN(lenA, lenB);
            }
            _TEMPLATE(T, poly_set_length) (G, lenG, ctx);

            if (G->length == 1)
                TEMPLATE(T, one) (G->coeffs, ctx);
            else
                TEMPLATE(T, poly_make_monic) (G, G, ctx);
        }
    }
}


#endif
