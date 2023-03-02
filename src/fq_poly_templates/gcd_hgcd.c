/*
    Copyright (C) 2011 William Hart
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

#include <stdlib.h>
#include "mpn_extras.h"

#define __set(B, lenB, A, lenA)                     \
do {                                                \
    _TEMPLATE(T, vec_set)((B), (A), (lenA), ctx);   \
    (lenB) = (lenA);                                \
} while (0)

#define __rem(R, lenR, A, lenA, B, lenB, invB)                              \
do {                                                                        \
    if ((lenA) >= (lenB))                                                   \
    {                                                                       \
        _TEMPLATE(T, poly_rem)((R), (A), (lenA), (B), (lenB), (invB), ctx); \
        (lenR) = (lenB) - 1;                                                \
        TEMPLATE(CAP_T, VEC_NORM)(R, lenR, ctx);                            \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        _TEMPLATE(T, vec_set)((R), (A), (lenA), ctx);                       \
        (lenR) = (lenA);                                                    \
    }                                                                       \
} while (0)

/*
    XXX: Incidentally, this implementation currently supports aliasing.  
    But since this may change in the future, no function other than 
    TEMPLATE(T, poly_gcd_hgcd)() should rely on this.
 */

slong
_TEMPLATE(T, poly_gcd_hgcd) (TEMPLATE(T, struct) * G,
                             const TEMPLATE(T, struct) * A, slong lenA,
                             const TEMPLATE(T, struct) * B, slong lenB,
                             const TEMPLATE(T, t) invB,
                             const TEMPLATE(T, ctx_t) ctx)
{
    slong cutoff, lenG, lenJ, lenR;
    TEMPLATE(T, struct) * J = _TEMPLATE(T, vec_init) (2 * lenB, ctx);
    TEMPLATE(T, struct) * R = J + lenB;
    TEMPLATE(T, t) inv;

    if (fmpz_bits(TEMPLATE(T, ctx_prime) (ctx)) <= 8)
        cutoff = TEMPLATE(CAP_T, POLY_SMALL_GCD_CUTOFF);
    else
        cutoff = TEMPLATE(CAP_T, POLY_GCD_CUTOFF);

    __rem(R, lenR, A, lenA, B, lenB, invB);

    if (lenR == 0)
    {
        __set(G, lenG, B, lenB);
    }
    else
    {
        TEMPLATE(T, init) (inv, ctx);

        _TEMPLATE(T, poly_hgcd) (NULL, NULL, G, &(lenG), J, &(lenJ), B, lenB,
                                 R, lenR, ctx);

        while (lenJ != 0)
        {
            TEMPLATE(T, inv) (inv, J + lenJ - 1, ctx);
            __rem(R, lenR, G, lenG, J, lenJ, inv);

            if (lenR == 0)
            {
                __set(G, lenG, J, lenJ);
                break;
            }
            if (lenJ < cutoff)
            {
                TEMPLATE(T, inv) (inv, R + lenR - 1, ctx);
                lenG =
                    _TEMPLATE(T, poly_gcd_euclidean) (G, J, lenJ, R, lenR, inv,
                                                      ctx);
                break;
            }

            _TEMPLATE(T, poly_hgcd) (NULL, NULL, G, &(lenG), J, &(lenJ), J,
                                     lenJ, R, lenR, ctx);
        }
        TEMPLATE(T, clear) (inv, ctx);
    }
    _TEMPLATE(T, vec_clear) (J, 2 * lenB, ctx);

    return lenG;
}

void TEMPLATE(T, poly_gcd_hgcd) (TEMPLATE(T, poly_t) G,
                                 const TEMPLATE(T, poly_t) A,
                                 const TEMPLATE(T, poly_t) B,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    if (A->length < B->length)
    {
        TEMPLATE(T, poly_gcd_hgcd) (G, B, A, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        TEMPLATE(T, poly_t) tG;
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
                TEMPLATE(T, poly_init2) (tG, FLINT_MIN(lenA, lenB), ctx);
                g = tG->coeffs;
            }
            else
            {
                TEMPLATE(T, poly_fit_length) (G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            TEMPLATE(T, init) (invB, ctx);
            TEMPLATE(T, inv) (invB, B->coeffs + lenB - 1, ctx);

            lenG = _TEMPLATE(T, poly_gcd_hgcd) (g, A->coeffs, lenA,
                                                B->coeffs, lenB, invB, ctx);

            if (G == A || G == B)
            {
                TEMPLATE(T, poly_swap) (tG, G, ctx);
                TEMPLATE(T, poly_clear) (tG, ctx);
            }
            G->length = lenG;

            if (G->length == 1)
                TEMPLATE(T, one) (G->coeffs, ctx);
            else
                TEMPLATE(T, poly_make_monic) (G, G, ctx);

            TEMPLATE(T, clear) (invB, ctx);
        }
    }
}

#undef __set
#undef __rem

#endif
