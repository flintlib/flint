/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

slong
_TEMPLATE(T, poly_gcd) (TEMPLATE(T, struct)*G,
                        const TEMPLATE(T, struct) * A, slong lenA,
                        const TEMPLATE(T, struct) * B, slong lenB,
                        const TEMPLATE(T, t) invB,
                        const TEMPLATE(T, ctx_t) ctx)
{
    slong cutoff;
    if (fmpz_bits(TEMPLATE(T, ctx_prime) (ctx)) <= 8)
        cutoff = TEMPLATE(CAP_T, POLY_SMALL_GCD_CUTOFF);
    else
        cutoff = TEMPLATE(CAP_T, POLY_GCD_CUTOFF);

    if (lenA < cutoff)
        return _TEMPLATE(T, poly_gcd_euclidean) (G, A, lenA, B, lenB, invB,
                                                 ctx);
    else
        return _TEMPLATE(T, poly_gcd_hgcd) (G, A, lenA, B, lenB, invB, ctx);
}

void
TEMPLATE(T, poly_gcd) (TEMPLATE(T, poly_t) G,
                       const TEMPLATE(T, poly_t) A,
                       const TEMPLATE(T, poly_t) B,
                       const TEMPLATE(T, ctx_t) ctx)
{
    if (A->length < B->length)
    {
        TEMPLATE(T, poly_gcd) (G, B, A, ctx);
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
            lenG = _TEMPLATE(T, poly_gcd) (g, A->coeffs, lenA,
                                           B->coeffs, lenB, invB, ctx);
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
