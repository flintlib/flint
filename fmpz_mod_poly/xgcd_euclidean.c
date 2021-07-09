/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

slong _fmpz_mod_poly_xgcd_euclidean(fmpz *G, fmpz *S, fmpz *T, 
                                   const fmpz *A, slong lenA, 
                                   const fmpz *B, slong lenB, 
                                   const fmpz_t invB, const fmpz_t p)
{
    _fmpz_vec_zero(G, lenB);
    _fmpz_vec_zero(S, lenB - 1);
    _fmpz_vec_zero(T, lenA - 1);

    if (lenB == 1)
    {
        fmpz_set(G + 0, B + 0);
        fmpz_one(T + 0);
        return 1;
    }
    else
    {
        fmpz *Q, *R;
        slong lenQ, lenR;
        TMP_INIT;
        TMP_START;

        FMPZ_VEC_TMP_INIT(Q, 2*lenA);
        R = Q + lenA;

        _fmpz_mod_poly_divrem(Q, R, A, lenA, B, lenB, invB, p);
        lenR = lenB - 1;
        FMPZ_VEC_NORM(R, lenR);

        if (lenR == 0)
        {
            _fmpz_vec_set(G, B, lenB);
            fmpz_one(T + 0);

            FMPZ_VEC_TMP_CLEAR(Q, 2*lenA);

            TMP_END;
            return lenB;
        }
        else
        {
            fmpz_t inv;
            fmpz *D, *U, *V1, *V3, *W;
            slong lenD, lenU, lenV1, lenV3, lenW;

            fmpz_init(inv);
            FMPZ_VEC_TMP_INIT(W, FLINT_MAX(5*lenB, lenA + lenB));
            D  = W  + lenB;
            U  = D  + lenB;
            V1 = U  + lenB;
            V3 = V1 + lenB;

            lenU = 0;
            _fmpz_vec_set(D, B, lenB);
            lenD = lenB;
            fmpz_one(V1 + 0);
            lenV1 = 1;
            lenV3 = 0;
            FMPZ_VEC_SWAP(V3, lenV3, R, lenR);

            do {
                fmpz_invmod(inv, V3 + (lenV3 - 1), p);
                _fmpz_mod_poly_divrem_basecase(Q, D, D, lenD, V3, lenV3, inv, p);
                lenQ = lenD - lenV3 + 1;
                lenD = lenV3 - 1;
                FMPZ_VEC_NORM(D, lenD);

                if (lenV1 >= lenQ)
                    _fmpz_mod_poly_mul(W, V1, lenV1, Q, lenQ, p);
                else
                    _fmpz_mod_poly_mul(W, Q, lenQ, V1, lenV1, p);
                lenW = lenQ + lenV1 - 1;

                _fmpz_mod_poly_sub(U, U, lenU, W, lenW, p);
                lenU = FLINT_MAX(lenU, lenW);
                FMPZ_VEC_NORM(U, lenU);

                FMPZ_VEC_SWAP(U, lenU, V1, lenV1);
                FMPZ_VEC_SWAP(D, lenD, V3, lenV3);

            } while (lenV3 != 0);

            _fmpz_vec_set(G, D, lenD);
            _fmpz_vec_set(S, U, lenU);
            {
                lenQ = lenA + lenU - 1;

                _fmpz_mod_poly_mul(Q, A, lenA, S, lenU, p);
                _fmpz_mod_poly_neg(Q, Q, lenQ, p);
                _fmpz_mod_poly_add(Q, G, lenD, Q, lenQ, p);

                _fmpz_mod_poly_divrem(T, W, Q, lenQ, B, lenB, invB, p);
            }

            FMPZ_VEC_TMP_CLEAR(W, FLINT_MAX(5*lenB, lenA + lenB));
            FMPZ_VEC_TMP_CLEAR(Q, 2*lenA);
            fmpz_clear(inv);

            TMP_END;
            return lenD;
        }
    }
}

void fmpz_mod_poly_xgcd_euclidean(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
                             fmpz_mod_poly_t T, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_xgcd_euclidean(G, T, S, B, A, ctx);
    }
    else  /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        const fmpz * p = fmpz_mod_ctx_modulus(ctx);
        fmpz_t inv;

        fmpz_init(inv);
        if (lenA == 0)  /* lenA = lenB = 0 */
        {
            fmpz_mod_poly_zero(G, ctx);
            fmpz_mod_poly_zero(S, ctx);
            fmpz_mod_poly_zero(T, ctx);
        }
        else if (lenB == 0)  /* lenA > lenB = 0 */
        {
            fmpz_invmod(inv, fmpz_mod_poly_lead(A, ctx), p);
            fmpz_mod_poly_scalar_mul_fmpz(G, A, inv, ctx);
            fmpz_mod_poly_zero(T, ctx);
            fmpz_mod_poly_set_fmpz(S, inv, ctx);
        }
        else  /* lenA >= lenB >= 1 */
        {
            fmpz *g, *s, *t;
            slong lenG;

            if (G == A || G == B)
            {
                g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }
            if (S == A || S == B)
            {
                s = _fmpz_vec_init(lenB);
            }
            else
            {
                fmpz_mod_poly_fit_length(S, lenB, ctx);
                s = S->coeffs;
            }
            if (T == A || T == B)
            {
                t = _fmpz_vec_init(lenA);
            }
            else
            {
                fmpz_mod_poly_fit_length(T, lenA, ctx);
                t = T->coeffs;
            }

            fmpz_invmod(inv, fmpz_mod_poly_lead(B, ctx), p);
            lenG = _fmpz_mod_poly_xgcd_euclidean(g, s, t,
                                     A->coeffs, lenA, B->coeffs, lenB, inv, p);

            if (G == A || G == B)
            {
                _fmpz_vec_clear(G->coeffs, G->alloc);
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
            }
            if (S == A || S == B)
            {
                _fmpz_vec_clear(S->coeffs, S->alloc);
                S->coeffs = s;
                S->alloc  = lenB;
            }
            if (T == A || T == B)
            {
                _fmpz_vec_clear(T->coeffs, T->alloc);
                T->coeffs = t;
                T->alloc  = lenA;
            }

            _fmpz_mod_poly_set_length(G, lenG);
            _fmpz_mod_poly_set_length(S, FLINT_MAX(lenB - lenG, 1));
            _fmpz_mod_poly_set_length(T, FLINT_MAX(lenA - lenG, 1));
            _fmpz_mod_poly_normalise(S);
            _fmpz_mod_poly_normalise(T);

            if (!fmpz_is_one(fmpz_mod_poly_lead(G, ctx)))
            {
                fmpz_invmod(inv, fmpz_mod_poly_lead(G, ctx), p);
                fmpz_mod_poly_scalar_mul_fmpz(G, G, inv, ctx);
                fmpz_mod_poly_scalar_mul_fmpz(S, S, inv, ctx);
                fmpz_mod_poly_scalar_mul_fmpz(T, T, inv, ctx);
            }
        }
        fmpz_clear(inv);
    }
}

