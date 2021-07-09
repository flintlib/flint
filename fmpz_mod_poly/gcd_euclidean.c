/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

slong _fmpz_mod_poly_gcd_euclidean(fmpz *G, const fmpz *A, slong lenA, 
                                           const fmpz *B, slong lenB, 
                                           const fmpz_t invB, const fmpz_t p)
{
    if (lenB == 1)
    {
        fmpz_one(G);
        return 1;
    }
    else  /* lenA >= lenB > 1 */
    {
        const slong lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
        fmpz_t invR3;
        fmpz *Q, *R1, *R2, *R3, *T, *W;
        slong lenR2, lenR3;
        TMP_INIT;

        TMP_START;

        FMPZ_VEC_TMP_INIT(W, lenW);
        Q  = W;
        R1 = W + FLINT_MAX(lenA - lenB + 1, lenB);
        R2 = R1 + lenA;
        R3 = R2 + lenB;

        _fmpz_mod_poly_divrem(Q, R1, A, lenA, B, lenB, invB, p);

        lenR3 = lenB - 1;
        FMPZ_VEC_NORM(R1, lenR3);

        if (lenR3 == 0)
        {
            _fmpz_vec_set(G, B, lenB);
            FMPZ_VEC_TMP_CLEAR(W, lenW);

            TMP_END;
            return lenB;
        }

        fmpz_init(invR3);

        T  = R3;
        R3 = R1;
        R1 = T;
        _fmpz_vec_set(R2, B, lenB);
        lenR2 = lenB;

        do
        {
            fmpz_invmod(invR3, R3 + (lenR3 - 1), p);

            _fmpz_mod_poly_divrem_basecase(Q, R2, R2, lenR2, R3, lenR3, invR3, p);
            lenR2 = lenR3 - 1;
            FMPZ_VEC_NORM(R2, lenR2);
            FMPZ_VEC_SWAP(R2, lenR2, R3, lenR3);
        } 
        while (lenR3 > 0);

        _fmpz_vec_set(G, R2, lenR2);

        FMPZ_VEC_TMP_CLEAR(W, lenW);
        fmpz_clear(invR3);

        TMP_END;
        return lenR2;
    }
}

void fmpz_mod_poly_gcd_euclidean(fmpz_mod_poly_t G, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_gcd_euclidean(G, B, A, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
        const fmpz * p = fmpz_mod_ctx_modulus(ctx);
        const slong lenA = A->length, lenB = B->length;
        slong lenG;
        fmpz *g;
    
        if (lenA == 0) /* lenA = lenB = 0 */
        {
            fmpz_mod_poly_zero(G, ctx);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            fmpz_mod_poly_make_monic(G, A, ctx);
        }
        else /* lenA >= lenB >= 1 */
        {
            fmpz_t invB;

            if (G == A || G == B)
            {
                g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            fmpz_init(invB);
            fmpz_invmod(invB, fmpz_mod_poly_lead(B, ctx), p);
            lenG = _fmpz_mod_poly_gcd_euclidean(g, A->coeffs, lenA,
                                                     B->coeffs, lenB, invB, p);
            fmpz_clear(invB);

            if (G == A || G == B)
            {
                _fmpz_vec_clear(G->coeffs, G->alloc);
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
                G->length = FLINT_MIN(lenA, lenB);
            }
            _fmpz_mod_poly_set_length(G, lenG);
    
            if (lenG == 1)
                fmpz_one(G->coeffs);
            else
                fmpz_mod_poly_make_monic(G, G, ctx);
        }
    }
}

