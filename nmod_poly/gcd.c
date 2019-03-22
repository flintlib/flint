/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"

slong _nmod_poly_gcd(mp_ptr G, mp_srcptr A, slong lenA, 
                              mp_srcptr B, slong lenB, nmod_t mod)
{
    const slong cutoff = FLINT_BIT_COUNT(mod.n) <= 8 ? 
                        NMOD_POLY_SMALL_GCD_CUTOFF : NMOD_POLY_GCD_CUTOFF;

    if (lenA < cutoff)
        return _nmod_poly_gcd_euclidean(G, A, lenA, B, lenB, mod);
    else
        return _nmod_poly_gcd_hgcd(G, A, lenA, B, lenB, mod);
}

void nmod_polydr_gcd(nmod_polydr_t G, 
            const nmod_polydr_t A, const nmod_polydr_t B, const nmod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        nmod_polydr_gcd(G, B, A, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_polydr_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_polydr_zero(G, ctx);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_polydr_make_monic(G, A, ctx);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_polydr_init2(tG, FLINT_MIN(lenA, lenB), ctx);
                g = tG->coeffs;
            }
            else
            {
                nmod_polydr_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, ctx->mod);

            if (G == A || G == B)
            {
                nmod_polydr_swap(tG, G, ctx);
                nmod_polydr_clear(tG, ctx);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_polydr_make_monic(G, G, ctx);
        }
    }
}

void nmod_poly_gcd(nmod_poly_t G, 
                             const nmod_poly_t A, const nmod_poly_t B)
{
    if (A->length < B->length)
    {
        nmod_poly_gcd(G, B, A);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        nmod_poly_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            nmod_poly_zero(G);
        } 
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            nmod_poly_make_monic(G, A);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                nmod_poly_init2(tG, A->mod.n, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                nmod_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd(g, A->coeffs, lenA,
                                               B->coeffs, lenB, A->mod);

            if (G == A || G == B)
            {
                nmod_poly_swap(tG, G);
                nmod_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                nmod_poly_make_monic(G, G);
        }
    }
}
