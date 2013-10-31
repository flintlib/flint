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
                       const fq_struct * B, long lenB, const fq_ctx_t ctx)
{
    long lenR1 = lenA, lenR2 = lenB, lenT;  /* assumes lenA \geq lenB */

    fq_struct *R1, *R2, *T;

    fq_t inv;

    fq_init(inv, ctx);
    T = _fq_vec_init(lenA, ctx);
    R1 = _fq_vec_init(lenA, ctx);
    R2 = _fq_vec_init(lenB, ctx);
    _fq_vec_set(R1, A, lenA, ctx);
    _fq_vec_set(R2, B, lenB, ctx);

    while (lenR2 > 0)
    {
        if (!fq_is_zero(R2 + (lenR2 - 1), ctx))
            fq_inv(inv, R2 + (lenR2 - 1), ctx);
        else
            flint_printf("something wrong with length");
        _fq_poly_rem(R1, R1, lenR1, R2, lenR2, inv, ctx);
        _fq_poly_normalise2(R1, &lenR1, ctx);

        _fq_vec_set(T, R1, lenR1, ctx);
        _fq_vec_set(R1, R2, lenR2, ctx);
        _fq_vec_set(R2, T, lenR1, ctx);
        lenT = lenR1;
        lenR1 = lenR2;
        lenR2 = lenT;
    }

    _fq_poly_set(G, R1, lenR1, ctx);
    _fq_vec_clear(T, lenA, ctx);
    _fq_vec_clear(R1, lenA, ctx);
    _fq_vec_clear(R2, lenB, ctx);
    fq_clear(inv, ctx);
    return lenR1;
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
        fq_poly_t tG;
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
                fq_poly_init2(tG, FLINT_MIN(lenA, lenB), ctx);
                g = tG->coeffs;
            }
            else
            {
                fq_poly_fit_length(G, FLINT_MIN(lenA, lenB), ctx);
                g = G->coeffs;
            }

            lenG = _fq_poly_gcd_euclidean(g, A->coeffs, lenA,
                                          B->coeffs, lenB, ctx);

            if (G == A || G == B)
            {
                fq_poly_swap(tG, G, ctx);
                fq_poly_clear(tG, ctx);
            }
            G->length = lenG;

            if (G->length == 1)
                fq_one(&(G->coeffs[0]), ctx);
            else
                fq_poly_make_monic(G, G, ctx);
        }
    }
}
