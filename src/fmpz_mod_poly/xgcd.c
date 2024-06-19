/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

slong _fmpz_mod_poly_xgcd(fmpz *G, fmpz *S, fmpz *T,
                                   const fmpz *A, slong lenA,
                                   const fmpz *B, slong lenB,
                                   const fmpz_t FLINT_UNUSED(invB), const fmpz_mod_ctx_t ctx)
{
    if (lenB == 1)
    {
        _fmpz_vec_zero(T, lenA - 1);
        fmpz_set(G + 0, B + 0);
        fmpz_one(T + 0);
        return 1;
    }
    else
    {
        slong lenG;
        gr_ctx_t gr_ctx;
        _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
                                           /* todo: use invB */

        if (FLINT_MIN(lenA, lenB) < FMPZ_MOD_POLY_GCD_CUTOFF)
            GR_MUST_SUCCEED(_gr_poly_xgcd_euclidean(&lenG, G, S, T, A, lenA, B, lenB, gr_ctx));
        else
            GR_MUST_SUCCEED(_gr_poly_xgcd_hgcd(&lenG, G, S, T, A, lenA, B, lenB, FMPZ_MOD_POLY_HGCD_CUTOFF, FMPZ_MOD_POLY_GCD_CUTOFF, gr_ctx));

        return lenG;
    }
}

void fmpz_mod_poly_xgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S,
                             fmpz_mod_poly_t T, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_xgcd(G, T, S, B, A, ctx);
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
            lenG = _fmpz_mod_poly_xgcd(g, s, t,
                                     A->coeffs, lenA, B->coeffs, lenB, inv, ctx);

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
