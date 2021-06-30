/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
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
_TEMPLATE(T, poly_xgcd_euclidean) (TEMPLATE(T, struct) * G,
                                   TEMPLATE(T, struct) * S,
                                   TEMPLATE(T, struct) * T,
                                   const TEMPLATE(T, struct) * A, slong lenA,
                                   const TEMPLATE(T, struct) * B, slong lenB,
                                   const TEMPLATE(T, t) invB,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE(T, vec_zero) (G, lenB, ctx);
    _TEMPLATE(T, vec_zero) (S, lenB - 1, ctx);
    _TEMPLATE(T, vec_zero) (T, lenA - 1, ctx);

    if (lenB == 1)
    {
        TEMPLATE(T, set) (G + 0, B + 0, ctx);
        TEMPLATE(T, one) (T + 0, ctx);
        return 1;
    }
    else
    {
        TEMPLATE(T, struct) * Q, *R;
        slong lenQ, lenR;

        Q = _TEMPLATE(T, vec_init) (2 * lenA, ctx);
        R = Q + lenA;

        _TEMPLATE(T, poly_divrem) (Q, R, A, lenA, B, lenB, invB, ctx);
        lenR = lenB - 1;
        TEMPLATE(CAP_T, VEC_NORM) (R, lenR, ctx);

        if (lenR == 0)
        {
            _TEMPLATE(T, vec_set) (G, B, lenB, ctx);
            TEMPLATE(T, one) (T + 0, ctx);

            _TEMPLATE(T, vec_clear) (Q, 2 * lenA, ctx);
            return lenB;
        }
        else
        {
            TEMPLATE(T, t) inv;
            TEMPLATE(T, struct) * D, *U, *V1, *V3, *W;
            slong lenD, lenU, lenV1, lenV3, lenW;

            TEMPLATE(T, init) (inv, ctx);
            W = _TEMPLATE(T, vec_init) (FLINT_MAX(5 * lenB, lenA + lenB), ctx);
            D = W + lenB;
            U = D + lenB;
            V1 = U + lenB;
            V3 = V1 + lenB;

            lenU = 0;
            _TEMPLATE(T, vec_set) (D, B, lenB, ctx);
            lenD = lenB;
            TEMPLATE(T, one) (V1 + 0, ctx);
            lenV1 = 1;
            lenV3 = 0;
            TEMPLATE(CAP_T, VEC_SWAP) (V3, lenV3, R, lenR);

            do
            {
                TEMPLATE(T, inv) (inv, V3 + (lenV3 - 1), ctx);
                _TEMPLATE(T, poly_divrem) (Q, R, D, lenD, V3, lenV3, inv, ctx);
                lenQ = lenD - lenV3 + 1;
                lenR = lenV3 - 1;
                TEMPLATE(CAP_T, VEC_NORM) (R, lenR, ctx);

                if (lenV1 >= lenQ)
                    _TEMPLATE(T, poly_mul) (W, V1, lenV1, Q, lenQ, ctx);
                else
                    _TEMPLATE(T, poly_mul) (W, Q, lenQ, V1, lenV1, ctx);
                lenW = lenQ + lenV1 - 1;

                _TEMPLATE(T, poly_sub) (U, U, lenU, W, lenW, ctx);
                lenU = FLINT_MAX(lenU, lenW);
                TEMPLATE(CAP_T, VEC_NORM) (U, lenU, ctx);

                TEMPLATE(CAP_T, VEC_SWAP) (U, lenU, V1, lenV1);
                {
                    TEMPLATE(T, struct) * __t;
                    slong __tn;

                    __t = D;
                    D = V3;
                    V3 = R;
                    R = __t;
                    __tn = lenD;
                    lenD = lenV3;
                    lenV3 = lenR;
                    lenR = __tn;
                }

            } while (lenV3 != 0);

            _TEMPLATE(T, vec_set) (G, D, lenD, ctx);
            _TEMPLATE(T, vec_set) (S, U, lenU, ctx);
            {
                lenQ = lenA + lenU - 1;

                _TEMPLATE(T, poly_mul) (Q, A, lenA, S, lenU, ctx);
                _TEMPLATE(T, poly_neg) (Q, Q, lenQ, ctx);
                _TEMPLATE(T, poly_add) (Q, G, lenD, Q, lenQ, ctx);

                _TEMPLATE(T, poly_divrem) (T, W, Q, lenQ, B, lenB, invB, ctx);
            }

            _TEMPLATE(T, vec_clear) (W, FLINT_MAX(5 * lenB, lenA + lenB), ctx);
            _TEMPLATE(T, vec_clear) (Q, 2 * lenA, ctx);
            TEMPLATE(T, clear) (inv, ctx);

            return lenD;
        }
    }
}

void
TEMPLATE(T, poly_xgcd_euclidean) (TEMPLATE(T, poly_t) G,
                                  TEMPLATE(T, poly_t) S, TEMPLATE(T, poly_t) T,
                                  const TEMPLATE(T, poly_t) A,
                                  const TEMPLATE(T, poly_t) B,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    if (A->length < B->length)
    {
        TEMPLATE(T, poly_xgcd_euclidean) (G, T, S, B, A, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        TEMPLATE(T, t) inv;

        TEMPLATE(T, init) (inv, ctx);
        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            TEMPLATE(T, poly_zero) (G, ctx);
            TEMPLATE(T, poly_zero) (S, ctx);
            TEMPLATE(T, poly_zero) (T, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            TEMPLATE(T, inv) (inv, TEMPLATE(T, poly_lead) (A, ctx), ctx);
            TEMPLATE3(T, poly_scalar_mul, T) (G, A, inv, ctx);
            TEMPLATE(T, poly_zero) (T, ctx);
            TEMPLATE3(T, poly_set, T) (S, inv, ctx);
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            TEMPLATE(T, inv)(inv, B->coeffs + 0, ctx);
            TEMPLATE3(T, poly_set, T) (T, inv, ctx);
            TEMPLATE(T, poly_one)(G, ctx);
            TEMPLATE(T, poly_zero)(S, ctx);
        }
        else                    /* lenA >= lenB >= 2 */
        {
            TEMPLATE(T, struct) * g, *s, *t;
            slong lenG;

            if (G == A || G == B)
            {
                g = _TEMPLATE(T, vec_init) (FLINT_MIN(lenA, lenB), ctx);
            }
            else
            {
                TEMPLATE(T, poly_fit_length) (G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }
            if (S == A || S == B)
            {
                s = _TEMPLATE(T, vec_init) (lenB, ctx);
            }
            else
            {
                TEMPLATE(T, poly_fit_length) (S, lenB, ctx);
                s = S->coeffs;
            }
            if (T == A || T == B)
            {
                t = _TEMPLATE(T, vec_init) (lenA, ctx);
            }
            else
            {
                TEMPLATE(T, poly_fit_length) (T, lenA, ctx);
                t = T->coeffs;
            }

            TEMPLATE(T, inv) (inv, TEMPLATE(T, poly_lead) (B, ctx), ctx);
            lenG = _TEMPLATE(T, poly_xgcd_euclidean) (g, s, t, A->coeffs, lenA,
                                                      B->coeffs, lenB, inv,
                                                      ctx);

            if (G == A || G == B)
            {
                _TEMPLATE(T, vec_clear) (G->coeffs, G->alloc, ctx);
                G->coeffs = g;
                G->alloc = FLINT_MIN(lenA, lenB);
		G->length = G->alloc;
            }
            if (S == A || S == B)
            {
                _TEMPLATE(T, vec_clear) (S->coeffs, S->alloc, ctx);
                S->coeffs = s;
                S->alloc = lenB;
		S->length = S->alloc;
            }
            if (T == A || T == B)
            {
                _TEMPLATE(T, vec_clear) (T->coeffs, T->alloc, ctx);
                T->coeffs = t;
                T->alloc = lenA;
		T->length = T->alloc;
            }

            _TEMPLATE(T, poly_set_length) (G, lenG, ctx);
            _TEMPLATE(T, poly_set_length) (S, FLINT_MAX(lenB - lenG, 1), ctx);
            _TEMPLATE(T, poly_set_length) (T, FLINT_MAX(lenA - lenG, 1), ctx);
            _TEMPLATE(T, poly_normalise) (S, ctx);
            _TEMPLATE(T, poly_normalise) (T, ctx);

            if (!TEMPLATE(T, is_one) (TEMPLATE(T, poly_lead) (G, ctx), ctx))
            {
                TEMPLATE(T, inv) (inv, TEMPLATE(T, poly_lead) (G, ctx), ctx);
                TEMPLATE3(T, poly_scalar_mul, T) (G, G, inv, ctx);
                TEMPLATE3(T, poly_scalar_mul, T) (S, S, inv, ctx);
                TEMPLATE3(T, poly_scalar_mul, T) (T, T, inv, ctx);
            }
        }
        TEMPLATE(T, clear) (inv, ctx);
    }
}

#endif
