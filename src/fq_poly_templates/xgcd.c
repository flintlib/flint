/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "gr_poly.h"

/* todo: the gr method don't accept the precomputed inverse */
slong
_TEMPLATE(T, poly_xgcd) (TEMPLATE(T, struct) * G,
                                   TEMPLATE(T, struct) * S,
                                   TEMPLATE(T, struct) * T,
                                   const TEMPLATE(T, struct) * A, slong lenA,
                                   const TEMPLATE(T, struct) * B, slong lenB,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    slong cutoff;
    slong lenG;

#if defined(FQ_NMOD_POLY_H) || defined(FQ_ZECH_POLY_H)
    if (FLINT_BIT_COUNT(TEMPLATE(T, ctx_prime)(ctx)) <= 8)
#else
    if (fmpz_bits(TEMPLATE(T, ctx_prime)(ctx)) <= 8)
#endif
        cutoff = TEMPLATE(CAP_T, POLY_SMALL_GCD_CUTOFF);
    else
        cutoff = TEMPLATE(CAP_T, POLY_GCD_CUTOFF);

    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    if (FLINT_MIN(lenA, lenB) < cutoff)
        GR_MUST_SUCCEED(_gr_poly_xgcd_euclidean(&lenG, G, S, T, A, lenA, B, lenB, gr_ctx));
    else
        GR_MUST_SUCCEED(_gr_poly_xgcd_hgcd(&lenG, G, S, T, A, lenA, B, lenB, TEMPLATE(CAP_T, POLY_HGCD_CUTOFF), cutoff, gr_ctx));

    return lenG;
}

void
TEMPLATE(T, poly_xgcd) (TEMPLATE(T, poly_t) G,
                                  TEMPLATE(T, poly_t) S, TEMPLATE(T, poly_t) T,
                                  const TEMPLATE(T, poly_t) A,
                                  const TEMPLATE(T, poly_t) B,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    if (A->length < B->length)
    {
        TEMPLATE(T, poly_xgcd) (G, T, S, B, A, ctx);
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
            lenG = _TEMPLATE(T, poly_xgcd) (g, s, t, A->coeffs, lenA,
                                                      B->coeffs, lenB,
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
