/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

slong _fmpz_mod_poly_gcd(fmpz *G, const fmpz *A, slong lenA,
                                  const fmpz *B, slong lenB,
                                  const fmpz_mod_ctx_t ctx)
{
    if (lenB == 1)
    {
        fmpz_one(G);
        return 1;
    }
    else  /* lenA >= lenB > 1 */
    {
        slong lenG;
        gr_ctx_t gr_ctx;
        _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);

        if (FLINT_MIN(lenA, lenB) < FMPZ_MOD_POLY_GCD_CUTOFF)
            GR_MUST_SUCCEED(_gr_poly_gcd_euclidean(G, &lenG, A, lenA, B, lenB, gr_ctx));
        else
            GR_MUST_SUCCEED(_gr_poly_gcd_hgcd(G, &lenG, A, lenA, B, lenB, FMPZ_MOD_POLY_HGCD_CUTOFF, FMPZ_MOD_POLY_GCD_CUTOFF, gr_ctx));

        return lenG;
    }
}

void fmpz_mod_poly_gcd(fmpz_mod_poly_t G, const fmpz_mod_poly_t A,
                             const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fmpz_mod_poly_gcd(G, B, A, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
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
            if (G == A || G == B)
            {
                g = _fmpz_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                fmpz_mod_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            lenG = _fmpz_mod_poly_gcd(g, A->coeffs, lenA, B->coeffs, lenB, ctx);

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
