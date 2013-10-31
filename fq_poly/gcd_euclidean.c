/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdlib.h>
#include "fq_poly.h"

long
_fq_poly_gcd_euclidean(fq_struct * G, const fq_struct * A, long lenA,
                       const fq_struct * B, long lenB, const fq_t invB,
                       const fq_ctx_t ctx)
{
    if (lenB == 1)
    {
        fq_one(G, ctx);
        return 1;
    }
    else  /* lenA >= lenB > 1 */
    {
        const slong lenW = FLINT_MAX(lenA - lenB + 1, lenB) + lenA + 2 * lenB;
        fq_t invR3;
        fq_struct *Q, *R1, *R2, *R3, *T, *W;
        slong lenR2, lenR3;

        W  = _fq_vec_init(lenW, ctx);
        Q  = W;
        R1 = W + FLINT_MAX(lenA - lenB + 1, lenB);
        R2 = R1 + lenA;
        R3 = R2 + lenB;

        _fq_poly_divrem(Q, R1, A, lenA, B, lenB, invB, ctx);

        lenR3 = lenB - 1;
        FQ_VEC_NORM(R1, lenR3, ctx);

        if (lenR3 == 0)
        {
            _fq_vec_set(G, B, lenB, ctx);
            _fq_vec_clear(W, lenW, ctx);
            return lenB;
        }

        fq_init(invR3, ctx);

        T  = R3;
        R3 = R1;
        R1 = T;
        _fq_vec_set(R2, B, lenB, ctx);
        lenR2 = lenB;

        do
        {
            fq_inv(invR3, R3 + (lenR3 - 1), ctx);

            _fq_poly_divrem(Q, R1, R2, lenR2, R3, lenR3, invR3, ctx);
            lenR2 = lenR3--;
            FQ_VEC_NORM(R1, lenR3, ctx);
            T = R2; R2 = R3; R3 = R1; R1 = T;
        } 
        while (lenR3 > 0);

        _fq_vec_set(G, R2, lenR2, ctx);

        _fq_vec_clear(W, lenW, ctx);
        fq_clear(invR3, ctx);

        return lenR2;
    }

}

void
fq_poly_gcd_euclidean(fq_poly_t G,
                      const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)
{
    if (A->length < B->length)
    {
        fq_poly_gcd_euclidean(G, B, A, ctx);
    }
    else                        /* lenA >= lenB >= 0 */
    {
        long lenA = A->length, lenB = B->length, lenG;
        fq_t invB;
        fq_struct *g;

        if (lenA == 0)          /* lenA = lenB = 0 */
        {
            fq_poly_zero(G, ctx);
        }
        else if (lenB == 0)     /* lenA > lenB = 0 */
        {
            fq_poly_make_monic(G, A, ctx);
        }
        else                    /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                g = _fq_vec_init(FLINT_MIN(lenA, lenB), ctx);
            }
            else
            {
                fq_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }
            
            fq_init(invB, ctx);
            fq_inv(invB, fq_poly_lead(B, ctx), ctx);
            lenG = _fq_poly_gcd_euclidean(g, A->coeffs, lenA,
                                          B->coeffs, lenB, invB, ctx);
            fq_clear(invB, ctx);

            if (G == A || G == B)
            {
                _fq_vec_clear(G->coeffs, G->alloc, ctx);
                G->coeffs = g;
                G->alloc = FLINT_MIN(lenA, lenB);
                G->length = FLINT_MIN(lenA, lenB);
            }
            _fq_poly_set_length(G, lenG, ctx);

            if (G->length == 1)
                fq_one(G->coeffs, ctx);
            else
                fq_poly_make_monic(G, G, ctx);
        }
    }
}
